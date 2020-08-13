from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner, TriContourSet, LinearTriInterpolator
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from osgeo import ogr
import geopandas as gpd
from shapely.geometry import Point, Polygon
from fiona.crs import from_epsg
from geopandas import GeoSeries


MA = 10  # 超参数，设置shape file分割阈值，增加阈值会降低分割性能，相反减小阈值提高分割性能，过低会导致无对象生成

shapefile = "slopeV.shp"
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(shapefile, 0)
layer = dataSource.GetLayer()

# 在下方设置选取固定范围阈值
layer.SetAttributeFilter("Velocity <= -110")
l = len(layer)
a = np.empty([l, 3])
i = 0
for feature in layer:
    a[i,0] = feature.GetField("X")
    a[i, 1] = feature.GetField("Y")
    a[i, 2] = feature.GetField("Velocity")
    # a[i, 2] = 1
    i = i+1
x_test=a[:,0]
y_test=a[:,1]
z_test=a[:,2]

tri = Triangulation(x_test, y_test)
# nedge = tri.edges.shape[0]
# l = tri.edges
#
# n = np.hypot(x_test[l[:,0]],x_test[l[:,1]])
# m = np.hypot(y_test[l[:,0]],y_test[l[:,1]])
# d = np.sqrt(n**2 + m**2)
a = tri.triangles
b1 = x_test[tri.triangles]
b2 = y_test[tri.triangles]
num = b1.shape[0]
c1 = np.empty([num,3])
c2 = np.empty([num,3])
for i in range (num):
    c1[i, 0] = b1[i, 0] - b1[i, 1]  # X
    c1[i, 1] = b1[i, 1] - b1[i, 2]
    c1[i, 2] = b1[i, 0] - b1[i, 2]

    c2[i, 0] = b2[i, 0] - b2[i, 1]  # Y
    c2[i, 1] = b2[i, 1] - b2[i, 2]
    c2[i, 2] = b2[i, 0] - b2[i, 2]

c1 = c1 / 1.67e-04
c2 = c2 / 1.67e-04
n =np.hypot(c1,c2).sum(axis=1)

# tri.set_mask(n>0.004)
tri.set_mask(n > MA)


fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
ax1.triplot(tri, 'bo-', lw=1)
tcf = ax1.tricontourf(tri, z_test)
fig1.colorbar(tcf)
ax1.tricontour(tri, z_test, colors='k')

ax1.set_title('triplot of Delaunay triangulation')
plt.show()


newdata = gpd.GeoDataFrame()
newdata['geometry'] = None
ntri = tri.triangles.shape[0]
X = tri.x
Y = tri.y
j = 0
for i in range(ntri):
    if (1-tri.mask[i]):
        n1 = tri.triangles[i,0]
        n2 = tri.triangles[i,1]
        n3 = tri.triangles[i,2]
        coordinates = [(X[n1], Y[n1]), (X[n2], Y[n2]), (X[n3], Y[n3])]
        # Create a Shapely polygon from the coordinate-tuple list
        poly = Polygon(coordinates)
        newdata.loc[j, 'geometry'] = poly
        j = j + 1

newdata = GeoSeries.set_crs(newdata, epsg=4326)  # 根据经纬度创建的矢量 需要定义WGS坐标系
newdata = newdata.to_crs("EPSG:4547")  # 为了计算面积需转换成投影坐标，这里转为CGCS 2000 114E 代号：4547

newdata['area'] = newdata.area
newdata.to_file("Senaatintori.shp")