import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
DTYPE = np.float64

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 14})

PI = np.pi

# mesh_file = np.genfromtxt("rect-0000-100x100"+".su2", DTYPE, skip_header=pow(100,2)+2, skip_footer=)

""" ----- Start ----- SU2 Mesh Input  """
# with open("SU2/mesh su2/rect-0000-100x100.su2", "r") as f:
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

# /* xOrigin, yOrigin, theta, lenght, width, Nx, Ny*/
# meshIdx = 0
meshNames = [# "singleBlock", 
            "bg", "comp"] 
mulFac = 5
mesh = np.array([(0.0, 0.0, 0.0 * PI/180.0, 1.0, 1.0, 10*mulFac, 10*mulFac),
        # (0.0, 0.0, 0.0 * PI/180.0, 1.0, 1.0, 100, 100),
        (0.4, 0.4, 1.0 * PI / 180.0, 0.4, 0.4, 4*mulFac, 4*mulFac)])

for meshIdx in range(2):
    meshDetails = mesh[meshIdx] 
    xOrgin, yOrigin = meshDetails[0], meshDetails[1]
    theta = meshDetails[2]
    L, W = meshDetails[3], meshDetails[4]
    Nx, Ny = int(meshDetails[5]), int(meshDetails[6])
    _x = np.linspace(0, L, Nx+2)
    _y = np.linspace(0, W, Ny+2)
    xLocal, yLocal = np.meshgrid(_x, _y)
    x = xOrgin + xLocal*np.cos(theta) - yLocal*np.sin(theta)
    y = yOrigin + xLocal*np.sin(theta) + yLocal*np.cos(theta)
        
    pointTypeData = np.loadtxt("Navier Stokes 2D/"+meshNames[meshIdx]+"PtType.txt")
    print(pointTypeData.shape)


    pointTypeLabels = ["unused", "calculated", "donor", "receiver", "donor buffer", "bc specified"];
    pointType = np.linspace(-0.5, 5.5, len(pointTypeLabels)+1)
    pointTypeLabels = [str(int(m))+" : "+n for m,n in zip(pointType+0.5, pointTypeLabels)]
    cmap = mpl.colormaps["turbo"]
    norm = colors.BoundaryNorm(pointType, cmap.N)


    imageScale = np.floor([Nx/mesh[1][-2], Ny/mesh[1][-1]])
    fig, ax = plt.subplots(figsize=(7 * imageScale[0], 5*imageScale[1]))
    im = ax.scatter(x, y, c=pointTypeData, cmap=cmap, norm=norm)
    # im = ax.plot(x, y)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_ticks(ticks=pointType[:-1]+0.5, labels=pointTypeLabels)
    # plt.show()
    fig.savefig(meshNames[meshIdx]+"HoleCut.png", bbox_inches ="tight", pad_inches = 0.2)

