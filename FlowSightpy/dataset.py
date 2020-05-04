
def load(datasets):
    '''Loads and merge all datasets

    datasets: a dictionary with names of datasts as heys, and file locations as values.

    return: pandas.DataFrame'''
    
    from functools import reduce
    import pandas as pd

    dflist = list()

    for name, file in datasets.items():
        df = pd.read_csv(file, sep="\t", skiprows=1, index_col=0)
        df["experiment"] = name
        df["old_index"] = df.index
        dflist.append(df)
    
    df = reduce(lambda x, y: pd.merge(x, y, how="outer"), dflist)

    return df

       
def to_adata(df, used_channels):
    '''Convert dataframe to an AnnData dataset with non flourescent channels as features. 
    All other channels are added as observations.

    Used_channels: describe which channels tif files have been exported for. These files should be stored under data/experiment/
    '''
    import numpy as np
    from sklearn.preprocessing import minmax_scale
    import scanpy as sc

    # Scale and normalize most features so they can be used with PCA
    df1 = df.copy()
    for column in df.columns: 
        
        if "Area" in column:
                df1[column] = np.log1p(df[column])
                
        if "Bright Detail Intensity" in column:
                df1[column] = np.log1p(df[column])
                
        if "Bkgd Mean" in column:
                df1[column] = minmax_scale(df[column])
                
        if "Contrast" in column:
                df1[column] = np.log1p(df[column])
                            
        if "Length" in column:
                df1[column] = np.log1p(df[column])
        
        if "Width" in column:
                df1[column] = np.log1p(df[column])
                
        if "Height" in column:
                df1[column] = np.log1p(df[column])
                
        if "Mean" in column:
                df1[column] = np.log1p(minmax_scale(df[column]))
                
        if "Median" in column:
                df1[column] = np.log1p(minmax_scale(df[column]))
    
    
    # save only non flourescent features (M01 and M06) to characterize the cells (keep the flourescent info for later)
    col_keep = [col for col in df1.columns if ("M01" in col or "M06" in col or "Ch01" in col or "Ch06" in col) and "Intensity" not in col]

    # Filter out all corrected features
    df_f = df1.filter(col_keep, axis=1)
    
    # Create adata with propper column and index + obs annotations
    adata = sc.AnnData(X=df_f.values)
    adata.var_names = df_f.columns
    adata.obs_names = df_f.index
    adata.obs = df1
    adata.obs["experiment"] = df["experiment"]
    adata.obs["old_index"] = df["old_index"]
    
    
    for ch in used_channels:
        adata.obs[ch] = [f"data/{e[1]['experiment']}/{e[1]['old_index']}_{ch}.ome.tif" for e in adata.obs.iterrows()]

    return adata


