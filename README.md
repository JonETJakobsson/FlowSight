# FlowSight
*Analyse data from the Amnis FlowSight to facilitate machine learning and 3rd party analysis.* 

Export features for each particle as .txt files using IDEAS software. Create a dictionary that links an "experiment name" to the exported .txt file. Create a folder called data/"experiment name"/ and export all particle images as tiff files to this directory for each experiment.
*Tip:* export only the channels that you need to save time.
This package uses [Scanpy](https://scanpy.readthedocs.io/en/latest/api/index.html) to perform dimentionality reduction and KNN based clustering, and [tifffile](https://github.com/cgohlke/tifffile) to read the ome.tif files.

Use `dataset.load()` to read in all data into a pandas dataframe. then `dataset.to_adata()` to convert this dataframe into a workable AnnData object. Only non flourescent features are used as variables. These are transformed and scaled to work with PCA. Other features are added as observations to facilitate explorative plotting. These observations are only transformed. Use `plot.PCA()` to run pca and plot an overview. Use `plot.cluster()` to find clusters in the data based on a neighberhood graph and the leiden algorithm. Use `plot.group_sample()` to get sample images of particles from specified groups. Use this to find out where your particles of interest are, and how many groups ofthese that you have.

See this "test/Test notebook.ipynb" for an example of all functions. 


