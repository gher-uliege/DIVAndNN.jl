# Read geotiff obtained from download_geotiff

import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal

imagefile = "foo3.tif"

gtif = gdal.Open(imagefile)

# info about the projection
print(gtif.GetProjectionRef())
arr = gtif.ReadAsArray()
trans = gtif.GetGeoTransform()
extent = (trans[0], trans[0] + gtif.RasterXSize*trans[1],
          trans[3] + gtif.RasterYSize*trans[5], trans[3])

# Example of plot
plt.imshow(arr, extent=extent,zorder=3,alpha=0.7)
plt.colorbar()
plt.show()

# Mask the "15" value
arrmasked = np.ma.masked_equal(arr, 15, copy=True)

# Plot with the masked field
plt.imshow(arrmasked, extent=extent, zorder=3, alpha=0.7)
plt.colorbar()
plt.show()
