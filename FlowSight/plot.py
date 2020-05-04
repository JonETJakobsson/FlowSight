

def PCA(adata, scale=True):
        '''Run filter, unit scaling and PCA to prepare for clustering'''
        import scanpy as sc
        #remove particles without masks
        sc.pp.filter_cells(adata, min_genes=2)

        # save the corrected feature values in raw for plotting
        adata.raw = adata.copy()
        # scale all features for dimentionality reduction and clustering
        if scale:
                sc.pp.scale(adata)
        
        sc.pp.pca(adata)
        sc.pl.pca_overview(adata, color="experiment")

        return adata


def cluster(adata, n_pcs=6, n_neighbors=10, resolution=1):

        import scanpy as sc
        # Create a KNN graph
        sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

        # Detect clusters using Leiden algorithm
        sc.tl.leiden(adata, resolution=resolution)

        # Embed KNN graph in 2 dimentions
        sc.tl.umap(adata)

        # Leiden groups and the 5 different experiments
        sc.pl.umap(adata, color=["leiden", "experiment"])

        return adata

def group_sample(adata, groupby="leiden", group="1", n_img=10, channels=["Ch1", "Ch6"]):
    '''Plot sample images of cells found in selected groups'''
    import numpy as np
    import matplotlib.pyplot as plt
    import tifffile
    
    adata_g = adata.copy()
    adata_g = adata_g[adata_g.obs[groupby] == group]
    df = adata_g.obs
    # get random samples from the group
    df = df.loc[np.random.choice(df.index, size=n_img, replace=False)]
    
    # Collect all images in a list
    im = list()
    for row in df.iterrows():
        files = list()
        for ch in channels:
                files.append(row[1][ch])
        im.append(tifffile.imread(files))


    #plot all images in rows
    i=1
    plt.figure(figsize=(len(channels)*2, n_img*1), facecolor='black')
    for row in im:
            for image in row:
                plt.subplot(len(im),len(channels), i, frameon=False)
                plt.imshow(image, cmap="gray")
                i += 1
    plt.subplots_adjust(hspace=0.2)

    