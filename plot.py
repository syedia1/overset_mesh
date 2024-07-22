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
    # fig = plt.figure(figsize=(6 + 3, 6 * l[0]/w[0] + 2))
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(6 + 3, 6 * 1 + 2))
    ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)
    for idx, mesh_name in enumerate(meshNames):
        meshDetails = mesh[idx]
        xOrigin, yOrigin, Theta, L, W, Nx, Ny = list(meshDetails)
        Nx, Ny = int(Nx)+2, int(Ny)+2
        print("Mesh '"+mesh_name+"' discretization- ", Nx, Ny, Theta * 180.0/PI, "deg")
        
        # Load data
        field_var = np.loadtxt("./SIMPLE/output_Upwind_"+mesh_name+"Mesh.dat", comments='#', dtype=DTYPE)
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
        pointMask = np.loadtxt("./SIMPLE/"+mesh_name+"PtType.txt")
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
        # # triang_Grid.set_mask(np.all(cellTypeMask[triang_Grid.triangles], axis=1)) # np.hypot(L/Nx, W/Ny)*1.1

        if idx == 0:
            var_max = np.sort(np.unique(field_var.flatten()))[-5]
            var_min = field_var.min()
            if k == 2:
                var_max = 0.8*var_max 
                # var_min = 0.8*var_min
                N_levels = N_levels * 5
            print(var_max, var_min)
        # Plotting field variable
        # k = 2 is pressure, which requries more contour lines for proper visulization

        LEVELS = np.linspace(var_min, var_max, N_levels)
        contour = plt.tricontour(triang_Grid, field_var.flatten(), levels = LEVELS, colors=colours[idx])
        plt.clabel(contour, inline = 1)
        
        # contour = plt.tricontour(triang_Grid, field_var.flatten(), levels = LEVELS)
        ax.add_artist(plt.Rectangle((xOrigin + L/Nx/2, yOrigin + W/Ny/2), L, W, angle = Theta * 180/PI,linewidth = 2, facecolor = "none", edgecolor=colours[idx]))
        # plt.xlabel(x_label)
        # plt.ylabel(y_label)
    # plt.colorbar(ticks = np.linspace(var_min, var_max, N_levels, endpoint=True), format = format)
    fig.suptitle(title + "-" + varNames[k])
    plt.savefig(filename+"_"+varNames[k])
    # if k == 2:
    #     plt.show()
    plt.close()

for k in range(3):
    fieldVarContour("NS contour for 2d plate", "lid_overset", k) 

def analytical_sol(xGrid, yGrid, Nx, Ny):
    T2, T1 = 100.0, 200.0
    L = l[0] + l[0]/nx[0] 
    W = w[0] + w[0]/ny[0]

    T_ana = np.zeros((Ny, Nx))
    for j in range(0, Ny):
        for i in range(0, Nx):
            x = xGrid[j][i]
            y = yGrid[j][i]
            fact = 0.0
            for n in range(1, 300):
                fact = fact + (pow(-1,n+1) + 1)/n * np.sin(n*PI*x/L) * np.sinh(n*PI*y/L) /np.sinh(n*PI*W/L)
            T_ana[j][i] = T1 + (T2-T1) * 2/PI * fact

    # return np.flip(T_ana)
    return T_ana

def error():

    fig, ax = plt.subplots(figsize=(6 + 3, 6 * l[0]/w[0] + 2))
    for idx, mesh_no in enumerate(["1", "2"]):
        # Load data
        Tn1 = np.loadtxt("num_mesh"+mesh_no+".txt", dtype=DTYPE) 
        L = l[idx]
        W = l[idx]
        Nx = nx[idx]
        Ny = ny[idx]

        # xAxis = np.linspace(0, L, Nx)+L/Nx/2 # dx/2 i.e. L/Nx/2 is added as we have values at cell center
        # yAxis = np.linspace(0, W, Ny)+W/Ny/2
        xAxis = np.linspace(0, L, Nx, dtype=DTYPE)
        yAxis = np.linspace(0, W, Ny, dtype=DTYPE)
        _xGrid, _yGrid = np.meshgrid(xAxis, yAxis)
        xGrid = x0[idx] + _xGrid * np.cos(theta[idx]) - _yGrid * np.sin(theta[idx])
        yGrid = y0[idx] + _xGrid * np.sin(theta[idx]) + _yGrid * np.cos(theta[idx])

        # T_analytical = analytical_sol(xGrid, yGrid, Nx, Ny)
        T_analytical = np.loadtxt("ana_mesh"+mesh_no+".txt", dtype=DTYPE) 

        error = np.linalg.norm(Tn1 - T_analytical, 2)
        print("Mesh "+mesh_no+" l2 error =  ", error)

        fig, ax = plt.subplots(figsize=(6 + 3, 6 * l[idx]/w[idx] + 2))
        LEVELS = np.linspace(var_min, var_max, N_levels, endpoint=True)
        
        contour1 = plt.tricontour(xGrid.flatten(), yGrid.flatten(), T_analytical.flatten(), levels = LEVELS, colors=colours[0])
        plt.clabel(contour1, inline = 1)
        ax.add_artist(plt.Rectangle((x0[idx], y0[idx]), L, W, angle = theta[idx]*180/PI,linewidth = 2, facecolor = "none", edgecolor=colours[1]))
        
        contour2 = plt.tricontour(xGrid.flatten(), yGrid.flatten(), Tn1.flatten(), levels = LEVELS, colors=colours[1])
        plt.clabel(contour2, inline = 1)
        # plt.show()
        plt.savefig("Temp analytical_"+mesh_no)
        plt.close()
    # plt.savefig("Temp analytical_")
    # plt.close()
        
        # plt.figure(figsize=(6 + 3, 6 * l[idx]/w[idx] + 2))
        # plt.scatter(xGrid, yGrid)
        # plt.savefig("scatter_"+mesh_no)
        # # plt.show()
        # plt.close()

# error()
