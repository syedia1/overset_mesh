import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

from scipy.interpolate import griddata

DTYPE = np.float64

PI = np.pi
# plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 14})

meshNames = [# "singleBlock", 
            "bg", "comp"] 
mulFac = 5
mesh = ((0.0, 0.0, 0.0 * PI/180.0, 1.0, 1.0, 10*mulFac, 10*mulFac),
        # (0.0, 0.0, 0.0 * PI/180.0, 1.0, 1.0, 100, 100),
        (0.4, 0.4, 0.0 * PI / 180.0, 0.4, 0.4, 4*mulFac, 4*mulFac))
varNames = ["U", "V", "P"]


def plot_Ghia(varName):
    if varName == "U":
        y = [1,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,0.5,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0]
        u = [1,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.2109,-0.15662,-0.1015,-0.06434,-0.04775,-0.04192,-0.03717,0]    
        axis = y
        var = u
    elif varName == "V":
        x = [0,0.061408882,0.072434916,0.094486983,0.15237366,0.221286371,0.24058193,0.499693721,0.802909648,0.858039816,0.902143951,0.940735069,0.954517611,0.965543645,1]
        v = [0.002061856,0.092783505,0.109278351,0.121649485,0.162886598,0.179381443,0.175257732,0.055670103,-0.245360825,-0.220618557,-0.167010309,-0.105154639,-0.084536082,-0.059793814,0.002061856]
        axis = x
        var = v
    plt.scatter(axis, var, c="black", marker="^", label="Ghia et al")

def centerline(title, filename, velComp):
    print("Potting - " + velComp)
    fig, ax = plt.subplots(figsize=(6 + 3, 6*1 + 2))
    # ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)  
    centerN = 100


    if velComp == "U":
        xCenterAxis = np.ones(1) * 0.5
        yCenterAxis = np.linspace(0, 1, centerN, endpoint=True)
    elif velComp == "V":
        xCenterAxis = np.linspace(0, 1, centerN, endpoint=True)
        yCenterAxis = np.ones(1) * 0.5
        

    xCenterGrid, yCenterGrid = np.meshgrid(xCenterAxis, yCenterAxis)

    k = varNames.index(velComp)
    plot_Ghia(varNames[k])

    for idx, mesh_name in enumerate(meshNames):
        meshDetails = mesh[idx]
        xOrigin, yOrigin, Theta, L, W, Nx, Ny = list(meshDetails)
        Nx, Ny = Nx+2, Ny+2
        print("Mesh '"+mesh_name+"' discretization- ", Nx, Ny, Theta * 180.0/PI, "deg")

        # Load data
        varData = np.loadtxt("./SIMPLE/output_Upwind_"+mesh_name+"Mesh.dat", comments='#', dtype=DTYPE)
        varData = varData[k*Ny:(k+1)*Ny, :]

        xAxis = np.linspace(0, L, Nx) 
        yAxis = np.linspace(0, W, Ny)
        _xGrid, _yGrid = np.meshgrid(xAxis, yAxis)
        _xGrid = _xGrid + L/Nx/2
        _yGrid = _yGrid + W/Ny/2
        xGrid = xOrigin + _xGrid * np.cos(Theta) - _yGrid * np.sin(Theta)
        yGrid = yOrigin + _xGrid * np.sin(Theta) + _yGrid * np.cos(Theta)

        pointMask = np.loadtxt("./SIMPLE/"+mesh_name+"PtType.txt")
        cellTypeMask = np.where((pointMask == 0), False, True)
        field_var = varData.flatten()
        field_var[np.invert(cellTypeMask)] = np.nan

        interpolatedCenterLine = griddata(np.stack((xGrid.flatten(), yGrid.flatten()), axis = 1), field_var, (xCenterGrid, yCenterGrid), method='linear')

        if velComp == "U":
            ax.plot(yCenterGrid.flatten(), interpolatedCenterLine.flatten(), label=mesh_name+" Mesh")
        elif velComp == "V":
            ax.plot(xCenterGrid.flatten(), interpolatedCenterLine.flatten(), label=mesh_name+" Mesh")

    plt.legend()
    fig.suptitle(title + " for " + velComp)
    ax.set_ylabel(velComp)
    if velComp == "U":
        ax.set_xlabel("y")
    elif velComp == "V":
        ax.set_xlabel("x")
        
    plt.savefig(filename+"_"+varNames[k])
    plt.close()

for var in varNames[:2]:
    centerline("Center line plot ", "lid_overset_centerline", var)