import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

DTYPE = np.float64

PI = np.pi

meshNames = ["bg", "comp"] 
mesh = np.array([(0.0, 0.0, 0.0, 1.0, 1.0, 50, 50), (0.445, 0.245, 45 * PI/180.0, 0.4, 0.4, 20, 20)])

l = [5, 1.5]
w = [5, 1.5]
mulFact = 2
nx = np.array([20, 6])*mulFact
ny = np.array([20, 6])*mulFact
x0 = [0, l[0]/3]
y0 = [0, w[0]/3]
theta = np.array([0, -25]) * PI/180.0 # counter-clock rotation
colours = ["k", "r", "b"]
var_max = 200
var_min = 100
format = "%0.1f"
N_levels = 20


print("Mesh discretization ", nx, ny, theta * 180.0/PI)
# # shift x0, y0 origins to cell center 
# for i in range(len(x0)):
#     dx = l[i]/nx[i]/2 # dx/2 i.e. L/Nx/2 is added as we have values at cell center
#     dy = w[i]/ny[i]/2
#     x0[i] = x0[i] + dx*np.cos(theta[i]) - dy*np.sin(theta[i])
#     y0[i] = y0[i] + dx*np.sin(theta[i]) + dy*np.cos(theta[i])

def temp_contour(title, filename, x_label = '$x$', y_label = '$y$'):
    # fig = plt.figure(figsize=(6 + 3, 6 * l[0]/w[0] + 2))
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(6 + 3, 6 * l[0]/w[0] + 2))
    ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)
    for idx, mesh_name in enumerate(meshNames):
        meshIdx = idx
        meshDetails = mesh[meshIdx]
        xOrgin, yOrigin = meshDetails[0], meshDetails[1]
        Theta = meshDetails[2]
        L, W = meshDetails[3], meshDetails[4]
        Nx, Ny = int(meshDetails[5]), int(meshDetails[6])
        x0[idx], y0[idx] = xOrgin, yOrigin
        theta[idx] = Theta

        # Load data
        Tn1 = np.loadtxt("num_mesh_"+mesh_name+".txt", dtype=DTYPE) 

        # L = l[idx]
        # W = l[idx]
        # Nx = nx[idx]
        # Ny = ny[idx]

        # Plotting temperature
        field_var = Tn1
        # field_var = np.array(field_var).T
        # field_var = np.flip(field_var, axis=0)
        xAxis = np.linspace(0, L, Nx) 
        yAxis = np.linspace(0, W, Ny)
        _xGrid, _yGrid = np.meshgrid(xAxis, yAxis)
        xGrid = x0[idx] + _xGrid * np.cos(theta[idx]) - _yGrid * np.sin(theta[idx])
        yGrid = y0[idx] + _xGrid * np.sin(theta[idx]) + _yGrid * np.cos(theta[idx])
        
        # remove unused cells which have -1 value
        xGrid, yGrid = xGrid.flatten(), yGrid.flatten()
        field_var = field_var.flatten()
        cellTypeMask = np.where(field_var < 0, False, True) 
        xGrid = xGrid[cellTypeMask]
        yGrid = yGrid[cellTypeMask]
        field_var = field_var[cellTypeMask]

        triang_Grid = tri.Triangulation(xGrid, yGrid)
        x_rem = xGrid[triang_Grid.triangles] - np.roll(xGrid[triang_Grid.triangles], 1, axis=1)
        y_rem = yGrid[triang_Grid.triangles] - np.roll(yGrid[triang_Grid.triangles], 1, axis=1)
        maxi = np.max(np.sqrt(x_rem**2 + y_rem**2), axis=1).reshape(-1)
        triang_Grid.set_mask((maxi > np.hypot(L/Nx, W/Ny)*1.1))



        # Plotting field variable
        LEVELS = np.linspace(var_min, var_max, N_levels, endpoint=True)
        # contour = plt.tricontour(xGrid.flatten(), yGrid.flatten(), field_var.flatten(), levels = LEVELS, colors=colours[idx])
        contour = plt.tricontour(triang_Grid, field_var.flatten(), levels = LEVELS, colors=colours[idx])
        plt.clabel(contour, inline = 1)
        ax.add_artist(plt.Rectangle((x0[idx], y0[idx]), L, W, angle = theta[idx]*180/PI,linewidth = 2, facecolor = "none", edgecolor=colours[idx]))
        
        # plt.xlabel(x_label)
        # plt.ylabel(y_label)
        # plt.show()
    # plt.colorbar(ticks = np.linspace(var_min, var_max, N_levels, endpoint=True), format = format)
    
    # ax.set_xlabel('100 C')
    # fig.supylabel('200 C', x=0.92, y=0.5)
    # fig.supxlabel('200 C', x=0.5, y=0.90)
    # ax.set_ylabel('200 C')
    # plt.show()
    plt.savefig(filename)
    fig.suptitle(title)
    plt.close()

temp_contour("Temp contour for 2d plate", "temp_overset") 
# # temp_contour(Tn1, "Temp contour for 2d plate mesh "+mesh_no, "temp_mesh"+mesh_no)

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
