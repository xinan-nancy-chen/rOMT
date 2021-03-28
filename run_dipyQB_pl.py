'''
import numpy as np
from nibabel import trackvis as tv
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle
from dipy.data import get_data
from dipy.viz import fvtk


fname = get_data('fornix')
streams, hdr = tv.read(fname)

streamlines = [i[0] for i in streams]
qb = QuickBundles(threshold=10.)
clusters = qb.cluster(streamlines)
print("Nb. clusters:", len(clusters))
print("Cluster sizes:", map(len, clusters))
print("Small clusters:", clusters < 10)
print("Streamlines indices of the first cluster:\n", clusters[0].indices)
print("Centroid of the last cluster:\n", clusters[-1].centroid)
ren = fvtk.ren()
ren.SetBackground(1, 1, 1)
fvtk.add(ren, fvtk.streamtube(streamlines, fvtk.colors.white))
fvtk.record(ren, n_frames=1, out_path='fornix_initial.png', size=(600, 600))

colormap = fvtk.create_colormap(np.arange(len(clusters)))

colormap_full = np.ones((len(streamlines), 3))
for cluster, color in zip(clusters, colormap):
    colormap_full[cluster.indices] = color

fvtk.clear(ren)
ren.SetBackground(1, 1, 1)
fvtk.add(ren, fvtk.streamtube(streamlines, colormap_full))
fvtk.record(ren, n_frames=1, out_path='fornix_clusters.png', size=(600, 600))
'''
import numpy as np 
#from nibabel import trackvis as tv
from dipy.segment.clustering import QuickBundles 
from dipy.io.pickles import save_pickle 
#from dipy.data import get_data
#from dipy.viz import fvtk
import scipy.io as spio 
from dipy.segment.metric import ResampleFeature 
from dipy.segment.metric import VectorOfEndpointsFeature 
from dipy.segment.metric import AveragePointwiseEuclideanMetric

mat = spio.loadmat('pl_cur.mat', squeeze_me=True) 
streams = mat['pl_cur'] 
feature = ResampleFeature(nb_points=124) 
metric = AveragePointwiseEuclideanMetric(feature=feature) 
qb = QuickBundles(threshold=4, metric=metric) 
pl = [i for i in streams] 
clusters = qb.cluster(pl) 
len(clusters) 
print("Nb. clusters:", len(clusters)) 
print("Cluster sizes:", list(map(len, clusters))) 
pli = [i.indices for i in clusters] 
pli_array = np.array([np.array(i) for i in pli]) 
spio.savemat('pli_array.mat', {'pli_array': pli_array}) 
pl_centroid = [i.centroid for i in clusters] 
pl_centroid_array = np.array([np.array(i) for i in pl_centroid]) 
spio.savemat('pl_centroid_array.mat', {'pl_centroid_array': pl_centroid_array})
