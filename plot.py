import matplotlib.pyplot as plt
import numpy as np

mesh_no = "2"
for mesh_no in ["1", "2"]:
    # Load data
    Tn1 = np.loadtxt("mesh"+mesh_no+".txt", dtype=float) 

    L = 1
    W = 1
    Nx, Ny = Tn1.shape
    Nx = Nx - 3
    Ny = Ny - 3
    var_max = 200
    var_min = 100

    # Plotting temperature
    domain_specs = [Nx, Ny, L, W]
    def temp_contour(domain_specs, field_var, title, filename, x_label = '$x$', y_label = '$y$'):
        field_var = np.array(field_var).T
        # field_var = np.flip(field_var, axis=0)
        Nx, Ny, L, W = domain_specs
        xGrid, yGrid = np.meshgrid(np.linspace(0, L, Nx), np.linspace(0, W, Ny))
        # print(xGrid.shape, yGrid.shape, field_var.shape)
        # Plotting field variable
        format = "%0.1f"

        plt.figure(figsize=(Nx//10 + 3, Ny//10 + 2))
        plt.tricontourf(xGrid.flatten(), yGrid.flatten(), field_var.flatten(), levels = np.linspace(var_min, var_max, 11, endpoint=True))
        plt.colorbar(ticks = np.linspace(var_min, var_max, 11, endpoint=True), format = format)
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        # plt.show()
        plt.savefig(filename)
        plt.close()

    temp_contour(domain_specs, Tn1[2:-1, 2:-1], "Temp contour for 2d plate mesh "+mesh_no, "temp_mesh"+mesh_no)