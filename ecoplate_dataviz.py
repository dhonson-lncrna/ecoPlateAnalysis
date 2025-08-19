import pandas as pd
import numpy as np

import os
import re
import sys
import yaml

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import math

def plot_curves(data,
                outdir,
                groupby,
                colorby,
                cmap='colorblind',):
    '''Creates scatter plots of blanked_590 values from Tidy.csv file.
    
    Params
    ______

    data : a Pandas DataFrame
        The data saved as [header]Tidy.csv from the data analysis script.

    outdir : str
        Directory where images will be saved.

    groupby : str or list
        Column(s) to plot in a single scatter. Each unique value will be 
        saved as an image file.

    colorby : str
        Column to use for coloring dots.

    cmap : str
        Colormap to use for coloring dots. Default: "colorblind"
    
    '''
    try:
        os.mkdir(outdir)
    except:
        pass
        
    if type(groupby) == str:
        ncol = groupby
    else:
        ncol = '_'.join(groupby)
        data[ncol] = data[groupby].agg(''.join, axis=1)

    for g in set(data[ncol]):
        subdf = data[data[ncol] == g]
        fig, ax  = plt.subplots(4,8, figsize=(24,12))

        for i, v in enumerate(set(subdf['metab'])):
            x = math.trunc(i / 8)
            y = int(i - 8*x)
            sns.scatterplot(data=subdf[subdf['metab'] == v],
                            x='hours',
                            y='blanked_590',
                            hue=colorby,
                            palette=cmap,
                            legend=False,
                            ax=ax[x,y])
            sns.regplot(data=subdf[subdf['metab'] == v],
                        x='hours',
                        y='blanked_590',
                        scatter=False,
                        color='black',
                        ci=90,
                        ax=ax[x,y])
            
            ax[x,y].set_ylabel(g+' blanked_590')
            ax[x,y].set_title(v)

        plt.tight_layout()
        plt.savefig(outdir+'/'+g+'_scatter.png',format='png',dpi=300)
        plt.close()
    return print(f'All scatters saved to {outdir}/')

def plot_pca(data,
             explained_variance,
             outdir,
             colorby,
             cmap='colorblind',
             show_cbar=False,
             show_plot=False):
    '''
    Function for plotting scatters of PCA data.

    Params
    ______

    data : a Pandas Dataframe
        Output of data analysis stored in [header]Pca.csv. The "sample"
        column should be split so that any value for "colorby" can be 
        recognized.

    explained_variance : list
        Explained variance by PC1 and PC2, output in [header]PcaReport.txt.
        Should be a list of two floats.

    outdir : str
        Directory where images will be saved.

    colorby : str or list
        Column(s) to use for coloring data. If list, the function will
        merge the columns and use the set of combined values.

    cmap : str
        Colormap to use for plotting.

    show_cbar : Bool
        If True, a colorbar will be added to the plot. Only useful for 
        non-categorical variables.

    show_plot : Bool
        Allows plot to be visualized in a Jupyter notebook if True. Default: False.
    '''
    try:
        os.mkdir(outdir)
    except:
        pass

    if type(colorby) == str:
        ncol = colorby
    else:
        ncol = '_'.join(colorby)
        data[ncol] = data[colorby].agg(''.join, axis=1)

    fig,ax = plt.subplots(1,1, figsize=(3.8,3))

    xlabel = f'PC1 ({explained_variance[0]})'
    ylabel = f'PC2 ({explained_variance[1]})'

    sns.scatterplot(data=data,
                    x='PC1',
                    y='PC2',
                    hue=ncol,
                    legend=True,
                    palette=cmap,
                    ax=ax)

    if show_cbar:
        ax.get_legend().remove()
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(ncol, rotation=270)
    else:
        ax.legend(bbox_to_anchor=(1,1))
        
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    plt.tight_layout()
    plt.savefig(''.join([outdir,'/pcaBy_',ncol,'.png']), format='png', dpi=300)
    if show_plot:
        plt.show()
    else:
        plt.close()

    return print(f'PCA plot saved to {outdir}/')
                

    
if __name__ == "__main__":
    # Check if a YAML file is provided as an argument
    if len(sys.argv) < 2:
        print("Please provide a YAML file path")
        sys.exit(1)
    
    # Get the YAML file path from command-line arguments
    yaml_file_path = sys.argv[1]

    with open(yaml_file_path, 'r') as file:
        params = yaml.safe_load(file)

    fdir = params['default_analysis']['fdir']

    # Grab files for analysis
    for f in os.listdir(fdir):
        if re.search(r'Tidy\.csv$',f):
            tidy = fdir+f
        elif re.search(r'Pca\.csv$',f):
            pca = fdir+f
        elif re.search(r'PcaReport\.txt$',f):
            pcareport = fdir+f
        else:
            pass

    # Import tidy data
    tidy_df = pd.read_csv(tidy)
    
    # Import PCA
    pca_df = pd.read_csv(pca)
    for i,v in enumerate(params['default_analysis']['sheet_info']):
        pca_df[v] = [j.split('_')[i] for j in pca_df['sample']]

    # Import explained variance
    with open(pcareport, 'r') as f:
        ls = f.readlines()

    evar= ls[1].split('[')[-1].split(']')[0].split(' ')
    explained_variance = []
    for i in evar:
        try:
            explained_variance.append(round(float(i), 3))
        except:
            pass

    # Generate scatters
    plot_curves(tidy_df,
                **params['scatters'])

    # Generate PCA
    plot_pca(pca_df,
             explained_variance,
             **params['pca'])
        