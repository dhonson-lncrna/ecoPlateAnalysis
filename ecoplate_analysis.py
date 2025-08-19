import pandas as pd
import numpy as np
import json
import os
import yaml
import sys

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def onefile_importer(fname,
                     sheet_json,
                     sample_info,
                     metab_json,
                     average_blanks=False,
                     zero_negatives=True,
                     plate_flip=False):
    '''
    Importer to convert EcoPlate data in Excel format to a tidy DataFrame.

    Params
    ______

    fname : str
        Path to an Excel file containing  EcoPlate data. The title of the files
        should be "[timepoint]h_[anyOtherInfo].xlsx". For example, 
        "24h_winter2025.xlsx". The function expects that every Excel file 
        represents a timepoint, and every sheet within the file represents a 
        plate. Plate data should be exported such that somewhere in the sheet is 
        a table with the columns "Well" (containing the well numbers) and "590" 
        (containing the absorbance values). It is not critical on which line the 
        table starts; the function will attempt to find a line that matches this
        format.

    sheet_json : str
        Path to a JSON file containing sample information for each sheet. The keys of 
        the JSON should be the sheet names (e.g., "Plate 1 - Sheet1") and the values
        should be the sample name (e.g., "25D_W_1"). If the sample name contains
        multiple pieces of relevant information (in this case condition, water/larvae,
        and replicate), the pieces should be separated by underscores. There are no limits
        on the amount of information you can store in the sample name, but all sample names
        must contain the same number of underscores.

    sample_info : tuple or None
        A tuple containing the desired column names for infomation derived from sample
        names. For samples with the format "25D_W_1", sample_info might be ("diapause",
        "water","replicate"). If None, the columns will be named "s1", "s2", etc.

    metab_json : str
        Path to a JSON file in which keys are the plate well positions and values are the 
        metabolites. The JSON provided in the "resources" folder of this repo is derived
        from this file: 
        https://www.biolog.com/wp-content/uploads/2023/08/00A-012-Rev-F-EcoPlate-IFU.pdf

    average_blanks : Bool
        If True, all Water well 590 values will be averaged and used as the blank for the 
        entire plate. If False, the A1 Water well will be used to blank columns 1-4, A5 
        for columns 5-8, and A9 for columns 9-12. Default: False
    
    zero_negatives : Bool
        If True, any blanked 590 absorbance values below zero will instead be set to zero/
        If False, negative values will remain negative. Default: True

    plate_flip : Bool
        Set to True if the plate was read such that H12 was accidentally put in the A1
        position. The function will remap all of the values to the correct wells before
        labeling metabolites.

    Output
    ______

    A Pandas DataFrame containing columns for well, metabolite, raw 590 value, blanked
    590 value, and any sample information.     
    
    '''
    # Get sample ID and metabolite info
    with open(sheet_json,'r') as f:
        pid = json.load(f)

    with open(metab_json,'r') as f:
        eco = json.load(f)

    # Identify blank wells
    water = [k for k,v in eco.items() if v == 'Water']

    # Get sheets
    sheets = pd.ExcelFile(fname).sheet_names

    # Extract data for each sheet
    df_ls = []
    for s in sheets:
        # Identify header row to get number of rows to skip
        df = pd.read_excel(fname,
                           sheet_name=s)
        skiprows = None
        for i in df.index:
            if 'Well' in list(df.loc[i,:]):
                skiprows=i+1
            else:
                pass

        # Verify that header was identified, otherwise raise an error
        try:
            int(skiprows)
        except:
            raise Exception('Header row not found for '+sheet_name)

        # Import tidy data
        df = pd.read_excel(fname,
                           sheet_name=s,
                           skiprows=skiprows)
        df = df.dropna(axis=1)
        df.columns = [str(i) for i in df.columns]

        # Correct plate flipping
        if plate_flip:
            orig_rows = [chr(65+i) for i in range(8)]
            orig_cols = [str(i) for i in np.arange(1,13)]
            
            df['orig_row'] = [i[0] for i in df['Well']]
            df['orig_col'] = [i[1:] for i in df['Well']]
            
            df['new_row'] = [orig_rows[::-1][orig_rows.index(r)] for r in df['orig_row']]
            df['new_col'] = [orig_cols[::-1][orig_cols.index(r)] for r in df['orig_col']]
            
            df['Well'] = df['new_row'] + df['new_col']
            df = df.drop(['orig_row','orig_col','new_row','new_col'], axis=1)
        else:
            pass 

        # Add metabolites
        df['metab'] = [eco[i] for i in df['Well']]

        # Add timepoints
        basename = fname.split('/')[-1]
        timepoint = basename.split('_')[0]
        df['timepoint'] = [timepoint for i in df.index]

        # Add sample information
        subdict = pid[timepoint]
        sample = subdict[s]
        subinfo = sample.split('_')

        df['sample'] = [sample for i in df.index]

        if sample_info:
            sdict = dict(zip(sample_info, subinfo))
        else:
            sdict = dict(zip(['s'+i for i,v in enumerate(subinfo)], 
                             subinfo))

        for col, info in sdict.items():
            df[col] = [info for i in df.index]

        #  Add blanks
        if average_blanks:
            avblank = np.mean(df[df['Well'].isin(water)]['590'])
            df['blanked_590'] = df['590'] - avblank
        else:
            df['column'] = [i[1:] for i in df['Well']]
            techreps = np.append(np.arange(1,12,4),13)
            for i, w in enumerate(water):
                blank = df.set_index('Well').loc[w,'590']
                cols = [str(i) for i in range(techreps[i],techreps[i+1])]
                subdf = df[df['column'].isin(cols)].copy()
                inds = subdf.index
                df.loc[inds,'blanked_590'] = subdf['590'] - blank
            df = df.drop('column', axis=1)
        # Zero negatives
        if zero_negatives:
            df['blanked_590'] = df['blanked_590'].clip(0)
        else:
            pass
                
        df_ls.append(df)

    # Concatenate, reindex
    df_all = pd.concat(df_ls)
    df_all.reset_index(drop=True, inplace=True)

    # Add hours column
    df_all['hours'] = [int(i[:-1]) for i in df_all['timepoint']]

    return df_all

def ecoplate_importer(fdir,
                      sheet_json,
                      sample_info,
                      metab_json,
                      average_blanks=False,
                      zero_negatives=True,
                      plate_flip=None):
    '''
    Importer to convert all EcoPlate data files in a folder to a tidy dataframe.

    Params
    ______

    fdir : str
        Path to the directory containing EcoPlate Excel files. The importer will attempt
        to import all Excel files in the folder. See "onefile_importer()" for additional
        details.

    sheet_json : str
        Path to a JSON file containing sample information for each sheet. The keys of 
        the JSON should be the sheet names (e.g., "Plate 1 - Sheet1") and the values
        should be the sample name (e.g., "25D_W_1"). If the sample name contains
        multiple pieces of relevant information (in this case condition, water/larvae,
        and replicate), the pieces should be separated by underscores. There are no limits
        on the amount of information you can store in the sample name, but all sample names
        must contain the same number of underscores.

    sample_info : tuple or None
        A tuple containing the desired column names for infomation derived from sample
        names. For samples with the format "25D_W_1", sample_info might be ("diapause",
        "water","replicate"). If None, the columns will be named "s1", "s2", etc.

    metab_json : str
        Path to a JSON file in which keys are the plate well positions and values are the 
        metabolites. The JSON provided in the "resources" folder of this repo is derived
        from this file: 
        https://www.biolog.com/wp-content/uploads/2023/08/00A-012-Rev-F-EcoPlate-IFU.pdf

    average_blanks : Bool
        If True, all Water well 590 values will be averaged and used as the blank for the 
        entire plate. If False, the A1 Water well will be used to blank columns 1-4, A5 
        for columns 5-8, and A9 for columns 9-12. Default: False
    
    zero_negatives : Bool
        If True, any blanked 590 absorbance values below zero will instead be set to zero/
        If False, negative values will remain negative. Default: True

    plate_flip : list or None
        A list of timepoints (in the format "24h") for which plates should be flipped.
        Function does not currently support plates with multiple orientations in the 
        same Excel file. Default: None.

    Output
    ______

    A Pandas DataFrame containing columns for well, metabolite, raw 590 value, blanked
    590 value, and any sample information for all files in the directory.
    
    '''
    # Extract files
    files = [i for i in os.listdir(fdir) if '.xlsx' in i]
    files = [i for i in files if '~' not in i]

    # Use onefile_importer() to get data for each file
    df_ls = []
    for f in files:
        if not plate_flip:
            pf = False
        elif f.split('_')[0] in plate_flip:
            pf = True
        else:
            pf = False
            
        df_ls.append(onefile_importer(fdir+f,
                                       sheet_json,
                                       sample_info,
                                       metab_json,
                                       average_blanks=False,
                                       zero_negatives=True,
                                       plate_flip=pf))

    return pd.concat(df_ls)

def averager(data,
             by,
             to_avg,
             keep = None,
             stdev = True,):
    '''
    A function to average imported EcoPlate data.

    Params
    ______

    data : a Pandas DataFrame
        The output of ecoplate_importer()

    by : str or list
        The columns that should be used to group entries for averaging.
        If a string, only that column will be used to group samples. If 
        a tuple or list, the union of the columns will be used for grouping.
        For example, if ['metab', 'sample'] is supplied for 'by', the values
        will be averaged by df['metab'] + df['sample']. An example value for
        this combination would be 'Water25EW1'.

    to_avg : str, tuple, or list
        Any values that should be averaged for the groupings identified by 
        'by'. Should be '590', 'blanked_590', or ['590', 'blanked_590'].

    keep : str, list, or None
        Any columns that are identical for all entries in the groups
        generated by 'by' and should be kept in the output dataframe.
        Usually this will be the 'hours' column, but could be others
        added to the dataframe after importing. Default: None.

    stdev : Bool
        Whether to make a column for standard deviation as well as 
        arithmetic mean. Default: True.
    '''

    # By can be one or more entries; handle that
    if type(by) != str:
        data[''.join(by)] = data[by].agg(''.join, axis=1)
        cols = by
        by = ''.join(by)
    else:
        cols = [by]

    # Keep is any information that is the same among all averaged samples, 
    # should not be included in the filtering to average, but should be retained
    # in the final output
    if keep:
        if type(keep) == str:
            cols.append(keep)
        else:
            cols += keep
    else:
        pass

    # Make a list to hold all the data. Iterate through all unique "by"
    # values and average desired values
    rows = []
    for b in set(data[by]):
        subdf = data[data[by] == b]
        # One or more values can be averaged 
        if type(to_avg) == str:
            avg = [np.mean(subdf[to_avg])]
            if stdev:
                avg.append(np.std(subdf[to_avg]))
            else:
                pass
        else:
            avg = [np.mean(subdf[a]) for a in to_avg]
            if stdev:
                avg += [np.mean(subdf[a]) for a in to_avg]
            else:
                pass
        
        out = [subdf.iloc[0][c] for c in cols] + avg
        rows.append(tuple(out))

    # Handle column names for averaging
    if type(to_avg) == str:
        cols.append(to_avg+'_mean')
        if stdev:
            cols.append(to_avg+'_std')
        else:
            pass
    else:
        cols += [a+'_mean' for a in to_avg]
        if stdev:
            cols += [a+'_std' for a in to_avg]
        else:
            pass
        
    return pd.DataFrame(rows, columns=cols) 

def integrator(av_df,
               by,
               xval,
               yval,
               exclude_negatives):
    '''
    Function to integrate average EcoPlate data over time. Calculates the area
    using the trapezoidal rule.

    Params
    ______

    av_df : a Pandas DataFrame
        The output of averager.

    by : str or list
        Which column(s) to group samples for integration. Each group
        should only have one average value for each timepoint.
        Typically should be the "by" from averager without "timepoint".

    xval : str
        Which column to use as x-values. Typically "hours".

    yval : str
        Which column to use as y-values. Typically "blanked_590_mean".

    exclude_negatives : str or None
        Parameter that instructs the function how to handle end blanked_590 values 
        that are lower than start blanked_590 values.

        flat : integrates by taking the t0 y-value and multiplying by the total 
            time.
        minimum : integrates by taking the minimum y-value and multiplying by the
            total time. Generally the safest option as it 
        drop : removes all values for metabolites where any end blanked_590 value 
            is higher than the start value.
        None : keeps all values and integrates as normal. Not recommended as this 
            will not distinguish between positive and negative trends during PCA.
        
    '''
    # Check validity of exclude_negatives
    if exclude_negatives in ['flat', 'minimum', 'drop', None]:
        pass
    else:
        raise ValueError('''exclude_negatives entry not valid. Choose from
                            ['flat', 'minimum', 'drop', None]''')
    # By can be one or more entries; handle that
    if type(by) != str:
        av_df['_'.join(by)] = av_df[by].agg('_'.join, axis=1)
        cols = by
        by = '_'.join(by)
    else:
        cols = [by]

    # Integrate
    rows = []
    skipped_samples = []
    for b in set(av_df[by]):
        subdf = av_df[av_df[by] == b].copy()
        subdf = subdf.sort_values(xval, ascending=True).reset_index(drop=True)
        # Exclude samples where endpoint lower than start point
        if subdf.loc[0,yval] > subdf.iloc[-1][yval]:
            skipped_samples.append(b)
            if exclude_negatives == 'flat':
                area = subdf.loc[0,yval] * (max(subdf[xval]) - min(subdf[xval]))
                out = [subdf.loc[0][c] for c in cols] + [area]
                rows.append(out)
            elif exclude_negatives == 'minimum':
                area = min(subdf[yval]) * (max(subdf[xval]) - min(subdf[xval]))
                out = [subdf.loc[0][c] for c in cols] + [area]
                rows.append(out)
            else:
                ind_areas = []
                for i in subdf.index[:-1]:
                    h = subdf.loc[i+1,xval] - subdf.loc[i,xval]
                    a = subdf.loc[i,yval]
                    b = subdf.loc[i+1,yval]
                    ind_areas.append(h * (a+b) / 2)
                area = np.sum(ind_areas)
                out = [subdf.loc[0][c] for c in cols] + [area]
                rows.append(out)     
        else:
            ind_areas = []
            for i in subdf.index[:-1]:
                h = subdf.loc[i+1,xval] - subdf.loc[i,xval]
                a = subdf.loc[i,yval]
                b = subdf.loc[i+1,yval]
                ind_areas.append(h * (a+b) / 2)
            area = np.sum(ind_areas)
            out = [subdf.loc[0][c] for c in cols] + [area]
            rows.append(out)

    cols.append('trapezoid_integration')
    outdf = pd.DataFrame(rows,columns=cols)
    if len(skipped_samples) == 0:
        pass
    elif exclude_negatives == 'drop':
        mind = by.split('_').index('metab')
        ex_metab = set([i.split('_')[mind] for i in skipped_samples])
        #print(len(ex_metab), ex_metab)
        outdf = outdf[~outdf['metab'].isin(ex_metab)]
    else:
        pass

    skipped_samples = [i+'\n' for i in skipped_samples]

    return outdf, skipped_samples

def metabolite_pca(int_df,
                   cols='metab',
                   index='sample',
                   n_components=2, 
                   scale_data=True):
    """
    Perform PCA on metabolite data and return results in long format.
    Written with assistance of Claude.ai.
    
    Params
    ______
    df : Pandas DataFrame 
        The output of integrator.

    cols : str
        Column in integrator dataframe to use as columns for PCA dataframe.
        Default: 'metab'.

    index : str
        Column in integrator dataframe to use as index for PCA dataframe.
        Default: 'sample'.
        
    n_components : int
        Number of principal components. Default: 2.
        
    scale_data : Bool 
        If True, data will be standardized using the formula
        z = (x - u) / s, where x is the original value, u is the
        mean of the data, and s is the standard deviation of the data.
        Essentially, this is required to make the data roughly Gaussian.
        Default: True.
    
    Returns:
    DataFrame with columns: metabolite, sample, PC1, PC2, ...
    """
    # Reshape dataframe and remove water (should be zero for all)
    keep = [cols,index,'trapezoid_integration']
    df = int_df.loc[:,keep]
    df = df.pivot(columns=cols,index=index,values='trapezoid_integration')
    df = df.drop('Water',axis=1)
    
    # Standardize the data
    if scale_data:
        scaler = StandardScaler()
        df_scaled = pd.DataFrame(
            scaler.fit_transform(df),
            index=df.index,
            columns=df.columns
        )
    else:
        df_scaled = df.copy()
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(df_scaled)
    
    # Create DataFrame with PC scores
    pc_columns = [f'PC{i+1}' for i in range(n_components)]
    pca_df = pd.DataFrame(
        pca_result,
        index=df.index,  # samples
        columns=pc_columns)
    
    # Get PCA loadings (how much each metabolite contributes to each PC)
    loadings = pd.DataFrame(
        pca.components_.T,
        index=df.columns,  # metabolites
        columns=pc_columns)

    pca_df = pca_df.reset_index(drop=False)
    
    return pca_df, pca, loadings.round(3)

def default_analysis(fdir,
                     sheet_json,
                     sheet_info,
                     metab_json,
                     plate_flip,
                     header='ecoPlate',
                     output='defaultOutput'):
    '''
    Runs above functions with most common parameters. Output of all
    functions are saved as csv files and input parameters are described
    in ecoplate_importer(). The only unique parameter is "output", the
    name of the directory in which to store the values.
    '''
    # Generate a folder to store output
    try:
        os.mkdir(output)
    except:
        pass

    # Full header
    full = '/'.join([output, header])
    # Import all data
    df = ecoplate_importer(fdir,
                           sheet_json,
                           sheet_info,
                           metab_json,
                           average_blanks=False,
                           zero_negatives=True,
                           plate_flip=plate_flip)
    df.to_csv(full+'Tidy.csv', index=False)

    # Average data
    av_df = averager(df,
                 ['metab','timepoint']+sheet_info,
                 'blanked_590',
                 keep=['sample','hours'])
    av_df.to_csv(full+'Average.csv',index=False)

    # Trapezoid integration
    int_df, skipped = integrator(av_df,
                                 ['metab','sample'],
                                 'hours',
                                 'blanked_590_mean',
                                  None)
    int_df.to_csv(full+'Integration.csv', index=False)
    with open(full+'SkippedSamples.txt', 'w') as f:
        f.writelines(skipped)

    # Principal component analysis
    pca_df, pca, loadings = metabolite_pca(int_df)
    pca_df.to_csv(full+'Pca.csv', index=False)
    loadings.to_csv(full+'Loadings.csv',index=True)
    
    with open(full+'PcaReport.txt','w') as f:
        f.write(f"PCA Results Summary:\n")
        f.write(f"Explained variance ratio: {pca.explained_variance_ratio_}\n")
        f.write(f"Total explained variance: {pca.explained_variance_ratio_.sum():.3f}\n")

    return print("Default pipeline executed without errors.")

if __name__ == "__main__":
    # Written with assistance of Claude.ai
    # Check if a YAML file is provided as an argument
    if len(sys.argv) < 2:
        print("Please provide a YAML file path")
        sys.exit(1)
    
    # Get the YAML file path from command-line arguments
    yaml_file_path = sys.argv[1]

    with open(yaml_file_path, 'r') as file:
        params = yaml.safe_load(file)

    # Run default pipeline if specified
    if params['shared_params']['useDefault']:
        default_analysis(**params['default_analysis'])

    # Otherwise run with individually specified parameters
    else:
        # Set up output parameters
        outdir = params['shared_params']['outdir']
        prefix = params['shared_params']['header']
        try:
            os.mkdir(outdir)
        except:
            pass
        header = '/'.join([outdir,prefix])

        # Import data
        df = ecoplate_importer(**params['ecoplate_importer'])
        df.to_csv(header+'Tidy.csv', index=False)

        # Average
        av_df = averager(df, **params['averager'])
        av_df.to_csv(header+'Average.csv', index=False)

        # Integrate
        int_df, skipped = integrator(av_df, **params['integrator'])
        int_df.to_csv(header+'Integration.csv', index=False)

        if len(skipped) > 0:
            with open(header+'SkippedSamples.txt', 'w') as f:
                f.writelines(skipped)
        else:
            pass

        # Principal component analysis
        result_df, pca, loadings = metabolite_pca(int_df, **params['metabolite_pca'])
        result_df.to_csv(header+'Pca.csv', index=False)
        loadings.to_csv(header+'Loadings.csv',index=True)
        
        with open(header+'PcaReport.txt','w') as f:
            f.write(f"PCA Results Summary:\n")
            f.write(f"Explained variance ratio: {pca.explained_variance_ratio_}\n")
            f.write(f"Total explained variance: {pca.explained_variance_ratio_.sum():.3f}\n")

        print('Full pipeline executed correctly.')