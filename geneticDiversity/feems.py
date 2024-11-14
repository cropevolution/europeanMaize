#Modified from Kerstin's script
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import statsmodels.api as sm
import cartopy.feature as cfeature
from sklearn.metrics.pairwise import haversine_distances
from scipy.spatial.distance import pdist, squareform
import pickle
from matplotlib.offsetbox import AnchoredText

# viz
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs

# feems
import feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz, Objective, query_node_attributes
from feems.cross_validation import run_cv, comp_mats

#for historical borders
import shapefile
from cartopy.feature import ShapelyFeature
from shapely.geometry import shape


# change matplotlib fonts
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.sans-serif"] = "Helvetica"


def cov_to_dist(S):
    s2 = np.diag(S).reshape(-1, 1)
    ones = np.ones((s2.shape[0], 1))
    D = s2 @ ones.T + ones @ s2.T - 2 * S
    return(D)

print(feems.__version__)

#set up datapath
data_path = "/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/feems.maizeEU/feems.maizeEU/FEEMS-maizeEU.all.filteredQual/maizeEU.all.filteredQual-FEEMS"
#print path
print(data_path)



# read the genotype data and mean impute missing data
(bim, fam, G) = read_plink("{}/Maize-0.8.pruned".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

# setup graph
data_path = "/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/feems.maizeEU/feems.maizeEU"
coord = np.loadtxt("{}/european.fixed.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/europeanMaize-2.outer".format(data_path))  # outer coordinates
grid_path = "{}/grid_100.shp".format(data_path)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)

# construct spatial graph object
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

#run cross validation
# define grids
# reverse the order of lambdas and alphas for warmstart
lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]

# run cross-validation
cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)

# average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])

print("This this the lamb with the best fit: " + str(lamb_cv))

##plot cv error acrosss runs
fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lamb_grid), mean_cv_err, "."); #plot lambda against the CV error
ax.set_xlabel("log10(lambda)"); #xlabel 
ax.set_ylabel("CV Error"); #ylabel
ax.axvline(np.log10(lamb_cv), color = "orange") #draw line at lowest CV point 
plt.show()
#save figure
plt.savefig('CV_error-0.8-outer2.pdf',bbox_inches='tight',pad_inches=0)

##
sp_graph.fit(lamb_cv)

###PLOT WITH HISTORICAL BORDERS
fig = plt.figure(dpi=300)
axx = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

#add the historical borders as a feature 
##Load your historical borders data using shapefile.Reader
shapefile_path = '/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/feems.maizeEU/feems.maizeEU/historical_borders_txt/provincemax.shp'
sf_reader = shapefile.Reader(shapefile_path)

##Extract the shapes and records from the shapefile
shapes = sf_reader.shapes()

##Convert the shapes to Shapely geometries
historical_geometries = [shape(s) for s in shapes]

##Convert the Shapely geometries to a ShapelyFeature
historical_feature = ShapelyFeature(historical_geometries, ccrs.PlateCarree(), facecolor='none', edgecolor='grey',linewidth=0.2)

# Use the Viz class from your_module
v = Viz(axx, sp_graph,
        projection=ccrs.PlateCarree(),  # Use PlateCarree projection for both the Viz class and subplot
        edge_width=0.5,
        edge_alpha=1,
        edge_zorder=100,
        sample_pt_size=5,
        sample_pt_color="black",
        obs_node_linewidth=0,
        obs_node_size=6,
        obs_node_color="silver",
        coastline_linewidth=0.2)

# Draw the basics
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)

#add additional features 
anchored_text = AnchoredText("λ with lowest error\nλ = {:.3f}\ncv error = {:.5f}".format(lamb_cv, mean_cv_err[np.argmin(mean_cv_err), 0]), prop=dict(fontsize=8,backgroundcolor='whitesmoke'),loc='lower right')
axx.add_artist(anchored_text) #add the textbox, use the dict function for format changes 
axx.add_feature(cfeature.LAKES,linewidth=0.5,color="powderblue") #add lakes to the map
axx.add_feature(cfeature.RIVERS,linewidth=0.5, color="powderblue") # add rivers to the map
#axx.add_feature(cfeature.BORDERS,linewidth=0.2, color="grey", linestyle='dotted') # add borders to the map 
axx.add_feature(cfeature.OCEAN, color="aliceblue") # add the ocean to the map 

# Add the historical borders to the map
axx.add_feature(historical_feature)

#add the colorbar
v.cbar_font_size = 5
v.cbar_orientation = "horizontal"
v.cbar_ticklabelsize = 5
v.cbar_loc= "lower center"
v.draw_edge_colorbar()

# Show the plot
plt.savefig('/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/feems.maizeEU/feems.maizeEU/feems_perfect-0.8-outer2-historicalv2.pdf',bbox_inches='tight', pad_inches=0 )
plt.savefig('/home/mtakou/Dropbox/margarita/posdocAGStetter/maizeEU/feems.maizeEU/feems.maizeEU/feems_perfect-0.8-outer2-historicalv2.jpeg',bbox_inches='tight', pad_inches=0 )