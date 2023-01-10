"""
Cahn-Hilliard equation with Neumann homogeneous conditions.

  phi'= gamma * Laplace(w)                          in the unit square
  w = - epsilon^2 * Laplace(phi) + F'(phi)          in the unit square
  grad(phi) * n = grad(w) * n = 0                   on the boundary
  phi = random data between -0.01 and 0.01          at t = 0

where F(phi) = (phi^2-1)^2.

We will comupute the energy functional

E = epsilon^2/2 * \int_\Omega |\nabla \phi|^2 + \int_\Omega (phi^2-1)^2

in each time step.

DG semidiscrete space scheme and Eyre semidicrete time scheme
"""

from fenics import *
import random
import numpy as np
import matplotlib.pyplot as plt


# >>>>>>>>>>> Global variables >>>>>>>>>>>>>>>
# Define the energy vector
E = []
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def refine_cell(cell,phi_0,mean_value,d=0.2):
    p = cell.midpoint()
    if abs(phi_0(p)-mean_value) < d:
        do_refine = True
    else:
        do_refine = False
    return do_refine


def CH_time_step(mesh, degree, t, dt, time_iter, phi_0):

    eps = 0.01
    gamma = 1.0
    sigma = Constant(4.0) # penalty parameter

    savepic = 0 # Indicates if pictures are saved or not

    print("dt = %f" %(dt))


    #  Internal parameters:
    #    TOL: error tolerance to stop the iteration.
    #    REFINE_RATIO: proportion of cells to refine.
    #    LEVEL_MAX: maximum number of iterations
    #
    TOL = 1.0e-04
    REFINE_RATIO = 0.25
    LEVEL_MAX = 5

    h_tol = 0.02
    h_min = h_tol
    level = 0

    #
    #  ADAPTIVE LOOP
    #
    while level<=LEVEL_MAX and h_min>=h_tol:

        if level>0:
            cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
            cell_markers.set_all(False)
            for cell in cells(mesh):
                cell_markers[cell] = refine_cell(cell,phi_0,mean_value=0.0,d=0.5)


            mesh = refine(mesh, cell_markers)

        h_min = mesh.hmin()
        print(type(h_min))
        print("h_min = %f (h_tol = %f)" %(h_min,h_tol))

        plot(mesh)
        plt.show()

        P = FiniteElement("DG", mesh.ufl_cell(), degree) # Space of polynomials
        W = FunctionSpace(mesh, MixedElement([P,P])) # Space of functions
        V = FunctionSpace(mesh, P)

        n = FacetNormal(mesh)
        h = CellDiameter(mesh)

        phi_0 = interpolate(phi_0, V)

        if time_iter==0 and level==0:
            c = plot(phi_0)
            plt.title("Condición inicial")
            plt.colorbar(c)
            plt.show()

        # Define variational problem
        u = TrialFunction(W) # Meaningless function used to define the variational formulation
        v = TestFunction(W) # Meaningless function used to define the variational formulation

        phi, w = split(u)
        barw, barphi = split(v)

        a1 = phi * barw * dx \
            + dt * gamma * (dot(grad(w),grad(barw)) * dx \
            - dot(avg(grad(w)),n('+'))*jump(barw) * dS \
            - dot(avg(grad(barw)),n('+'))*jump(w) * dS \
            + sigma/h('+') * dot(jump(w), jump(barw)) * dS)
        L1 = phi_0 * barw * dx

        a2 = w * barphi * dx \
            - pow(eps,2) * (dot(grad(phi),grad(barphi))*dx \
            - dot(avg(grad(phi)),n('+'))*jump(barphi) * dS \
            - dot(avg(grad(barphi)),n('+'))*jump(phi) * dS \
            + sigma/h('+') * dot(jump(phi), jump(barphi)) * dS) \
            - 2 * phi * barphi * dx
        L2 = pow(phi_0,3) * barphi * dx - 3 * phi_0 * barphi * dx

        a = a1 + a2
        L = L1 + L2

        u = Function(W)

        # Compute solution
        solve(a == L, u)

        phi, w = u.split(True)

        level += 1

        # if h_min<h_tol:
        #     print("Se alcanzó la tolerancia")
        #     break

    # Update previous solution
    phi_0.assign(phi)

    # <<< End refinemnent loop

    c = plot(phi_0)
    plt.title("Solución en la etapa de refinado %d" % (level))
    plt.colorbar(c)
    plt.show()

    print('max = %f' % (phi_0.vector().get_local().max()))
    print('min = %f' % (phi_0.vector().get_local().min()))
    print('mass = %f' % (assemble(phi_0*dx)))


    if time_iter == 0:
        energy = assemble(0.5*pow(eps,2)*(dot(grad(phi_0),grad(phi_0))*dx - 2.0 * dot(avg(grad(phi_0)),n('+'))*jump(phi_0) * dS  + sigma/h('+') * pow(jump(phi_0),2) * dS) + 0.25 * pow(pow(phi_0,2)-1,2)*dx)
        E.append(energy)
        print('E =',energy)

    # Compute the energy
    energy = assemble(0.5*pow(eps,2)*(dot(grad(phi),grad(phi))*dx - 2.0 * dot(avg(grad(phi)),n('+'))*jump(phi) * dS  + sigma/h('+') * pow(jump(phi),2) * dS) + 0.25 * pow(pow(phi,2)-1,2)*dx)
    E.append(energy)
    print('E =',energy)

    if(savepic):
        if(i==0 or i==(num_steps/2-1) or i==(num_steps-1)):
            pic = plot(phi)
            plt.title("Función de campo de fase en t = %.4f" %(t))
            plt.colorbar(pic)
            plt.savefig("fig/DG-Eyre_nt-%d_t-%.4f.png" %(num_steps,t))
            plt.close()

    # Compute the mass
    print('mass = %f' % (assemble(phi*dx)))

    # Return the solution
    return phi

# "Main function"

if ( __name__ == '__main__' ):

    # Create mesh and define function space
    nx = ny = 32 # Boundary points
    print("nx = ny = %d" %(nx))
    mesh = UnitSquareMesh(nx,ny)

    # Random initial data
    degree = 1
    phi_0 = Expression("pow(x[0]-0.5,2) + pow(x[1]-0.5,2)<0.04? 1 : -1", degree=degree)

    T = 0.05            # final time
    num_steps = 100     # number of time steps
    dt = T / num_steps # time step size
    t = 0

    print("Iteraciones:")

    for i in range(num_steps):

        print("\nIteración %d:" %(i+1))

        # Update current time
        t += dt
        phi_0 = CH_time_step(mesh, degree, t, dt, i, phi_0)

    plt.plot(np.linspace(0,T,num_steps+1),E, color='red')
    plt.title("Energía discreta")
    plt.xlabel("Tiempo")
    plt.ylabel("Energía")
    if(savepic):
        plt.savefig("fig/DG-Eyre_nt-%d_energia.png" %(num_steps))
    plt.show()
