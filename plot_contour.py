import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

DTYPE = np.float64

PI = np.pi
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 14})

meshNames = [#"singleBlock", 
            "bg", "comp"] 
mulFac = 5
mesh = ((0.0, 0.0, 0.0 * PI/180.0, 1.0, 1.0, 10*mulFac, 10*mulFac),
        # (0.0, 0.0, 0.0 * PI/180.0, 1.0, 1.0, 100, 100),
        (0.42, 0.42, 0.0 * PI / 180.0, 0.4, 0.4, 6*mulFac, 6*mulFac))
colours = ["k", "r", "b"]
format = "%0.1f"
N_levels = 20
varNames = ["U", "V", "P"]

def fieldVarContour(title, filename, k=0):
    print("Potting - " + varNames[k])
    global N_levels

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(6 + 3, 6 * 1 + 2))
    ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)
    for idx, mesh_name in enumerate(meshNames):
        meshDetails = mesh[idx]
        xOrigin, yOrigin, Theta, L, W, Nx, Ny = list(meshDetails)
        Nx, Ny = int(Nx)+2, int(Ny)+2
        print("Mesh '"+mesh_name+"' discretization- ", Nx, Ny, Theta * 180.0/PI, "deg")
        
        # Load data
        field_var = np.loadtxt("Navier Stokes 2D/output_Upwind_"+mesh_name+"Mesh.dat", comments='#', dtype=DTYPE)
        # k = 0
        field_var = field_var[k*Ny:(k+1)*Ny, :]
        # print(Tn1.shape) 

        # generate the mesh for plotting
        xAxis = np.linspace(0, L, Nx) 
        yAxis = np.linspace(0, W, Ny)
        _xGrid, _yGrid = np.meshgrid(xAxis, yAxis)
        _xGrid = _xGrid + L/Nx/2
        _yGrid = _yGrid + W/Ny/2
        xGrid = xOrigin + _xGrid * np.cos(Theta) - _yGrid * np.sin(Theta)
        yGrid = yOrigin + _xGrid * np.sin(Theta) + _yGrid * np.cos(Theta)
        
        # remove masked/unused cells
        pointMask = np.loadtxt("Navier Stokes 2D/"+mesh_name+"PtType.txt")
        xGrid, yGrid = xGrid.flatten(), yGrid.flatten()
        field_var = field_var.flatten()
        cellTypeMask = np.where((pointMask == 0), False, True) 
        xGrid = xGrid[cellTypeMask]
        yGrid = yGrid[cellTypeMask]
        field_var = field_var[cellTypeMask]

        triang_Grid = tri.Triangulation(xGrid, yGrid)
        x_rem = xGrid[triang_Grid.triangles] - np.roll(xGrid[triang_Grid.triangles], 1, axis=1)
        y_rem = yGrid[triang_Grid.triangles] - np.roll(yGrid[triang_Grid.triangles], 1, axis=1)
        maxi = np.max(np.sqrt(x_rem**2 + y_rem**2), axis=1).reshape(-1)
        triang_Grid.set_mask((maxi > np.hypot(L/Nx, W/Ny)*1.1))

        if idx == 0:
            var_max = np.sort(np.unique(field_var.flatten()))[-5]
            var_min = field_var.min()
            # k = 2 is pressure, which requries more contour lines for proper visulization
            if k == 2:
                var_max = 0.8*var_max 
                N_levels = N_levels * 5
            print(var_max, var_min)
        # Plotting field variable

        LEVELS = np.linspace(var_min, var_max, N_levels)
        contour = plt.tricontour(triang_Grid, field_var.flatten(), levels = LEVELS, colors=colours[idx])
        plt.clabel(contour, inline = 1)
        

        ax.add_artist(plt.Rectangle((xOrigin + L/Nx/2, yOrigin + W/Ny/2), L, W, angle = Theta * 180/PI,linewidth = 2, facecolor = "none", edgecolor=colours[idx]))
    fig.suptitle(title + "-" + varNames[k])
    plt.savefig(filename+"_"+varNames[k])
    plt.close()

for k in range(3):
    fieldVarContour("NS contour for 2d plate", "lid_overset", k) 