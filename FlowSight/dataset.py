from functools import reduce
from sklearn.preprocessing import minmax_scale
import numpy as np

def load(dataset):
    '''Loads and merge all datasets

    datasets: a dictionary with names of datasts as heys, and file locations as values.

    return: pandas.DataFrame'''

    
    dflist = list()

    for name, file in datasets.items():
        df = pd.read_csv(file, sep="\t", skiprows=1, index_col=0)
        df["experiment"] = name
        df["old_index"] = df.index
        dflist.append(df)
    
    df = reduce(lambda x, y: pd.merge(x, y, how="outer"), dflist)

    return df

def scale(df):
    '''Scale the features of the dataset depending on type of measurement
    
    Fix scales and logs of the features (trying to approximate a normal distribution)
    '''
    
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
    return df1       
            

def to_adata(df, use_channels):
    # save only non flourescent features (M01 and M06) to characterize the cells (keep the flourescent info for later)
    col_keep = [col for col in df.columns if ("M01" in col or "M06" in col or "Ch01" in col or "Ch06" in col) and "Intensity" not in col]
    # Filter out all corrected features
    df_f = df1.filter(col_keep, axis=1,)
    
    # Create adata with propper column and index + obs annotations
    adata = sc.AnnData(X=df_f.values)
    adata.var_names = df_f.columns
    adata.obs_names = df_f.index
    adata.obs["tdTomato"] = df1["Mean Pixel_M04_Ch04"]
    adata.obs["experiment"] = df1["experiment"]
    adata.obs["old_index"] = df1["old_index"]

    #Set channel paths for image aquisition
    channels = ["Ch1", "Ch4", "Ch6"]
    for ch in channels:
        adata.obs[ch] = [f"data/{e[1][1]}/{e[1][2]}_{ch}.ome.tif" for e in adata.obs.iterrows()]
