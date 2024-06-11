from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
DTYPE = np.float64

# mesh_file = np.genfromtxt("rect-0000-100x100"+".su2", DTYPE, skip_header=pow(100,2)+2, skip_footer=)

""" ----- Start ----- SU2 Mesh Input  """
# with open("/Users/zlatangg/Documents/Overset/overset_mesh/SU2/mesh su2/rect-0000-100x100.su2", "r") as f:
#   nDim, nElem, nPoint = 0,0,0
#   lines = f.readlines()
#   for i in range(len(lines)):
#     if lines[i].startswith('%'): continue
#     elif lines[i].startswith('NDIME='): 
#       nDim = int(lines[i].split()[-1])
#     elif lines[i].startswith("NELEM="):
#       nElem = int(lines[i].split()[-1])
#       elemConn = np.genfromtxt(lines[i+1:i+1+nElem])[:, 1:-1]
#     elif lines[i].startswith("NPOIN="):
#       nPoint = int(lines[i].split()[-1])
#       meshPoint = np.genfromtxt(lines[i+1:i+1+nPoint])
#     elif lines[i].startswith("NMARK="):
#       continue
#     else:
#       continue
    
# print(meshPoint.shape)
# print(elemConn.shape)
# x = meshPoint[:, 0]
# y = meshPoint[:, 1]
""" ----- End ----- SU2 Mesh Input  """

nx, ny = (100 + 1, 100 + 1)
_x = np.linspace(0, 1, nx)
_y = np.linspace(0, 1, ny)
x, y = np.meshgrid(_x, _y)

pointTypeData = np.loadtxt("/Users/zlatangg/Documents/Overset/overset_mesh/bgPtType.txt")
print(pointTypeData.shape)


pointType = np.linspace(-0.5, 4.5, 6);
pointTypeLabels = ["unused", "calculated", "donor", "receiver", "donor buffer"];
pointTypeLabels = [str(int(m))+" : "+n for m,n in zip(pointType+0.5, pointTypeLabels)]
cmap = plt.cm.jet
norm = colors.BoundaryNorm(pointType, cmap.N)


fig, ax = plt.subplots()
im = ax.scatter(x, y, c=pointTypeData, cmap=cmap, norm=norm)
cbar = fig.colorbar(im, ax=ax)
cbar.set_ticks(ticks=pointType[:-1]+0.5, labels=pointTypeLabels)


plt.show()

