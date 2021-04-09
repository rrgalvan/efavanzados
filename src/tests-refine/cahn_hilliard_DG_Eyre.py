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

    #
    #  ADAPTIVE LOOP
    #
    for level in range(LEVEL_MAX):

        plot(mesh)
        plt.show()

        print("h = %f" %(mesh.hmax()))

        P = FiniteElement("DG", mesh.ufl_cell(), degree) # Space of polynomials
        W = FunctionSpace(mesh, MixedElement([P,P])) # Space of functions
        V = FunctionSpace(mesh, P)

        n = FacetNormal(mesh)
        h = CellDiameter(mesh)

        phi_n = interpolate(phi_0, V)

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
        L1 = phi_n * barw * dx

        a2 = w * barphi * dx \
            - pow(eps,2) * (dot(grad(phi),grad(barphi))*dx \
            - dot(avg(grad(phi)),n('+'))*jump(barphi) * dS \
            - dot(avg(grad(barphi)),n('+'))*jump(phi) * dS \
            + sigma/h('+') * dot(jump(phi), jump(barphi)) * dS) \
            - 2 * phi * barphi * dx
        L2 = pow(phi_n,3) * barphi * dx - 3 * phi_n * barphi * dx

        a = a1 + a2
        L = L1 + L2

        u = Function(W)

        # Compute solution
        solve(a == L, u)

        phi, w = u.split(True)

    # Update previous solution
    phi_n.assign(phi)

    # <<< End refinemnent loop

    c = plot(phi_n)
    plt.title("Condición inicial")
    plt.colorbar(c)
    plt.show()

    print('max = %f' % (phi_n.vector().get_local().max()))
    print('min = %f' % (phi_n.vector().get_local().min()))
    print('mass = %f' % (assemble(phi_n*dx)))


    if time_iter == 0:
        energy = assemble(0.5*pow(eps,2)*(dot(grad(phi_n),grad(phi_n))*dx - 2.0 * dot(avg(grad(phi_n)),n('+'))*jump(phi_n) * dS  + sigma/h('+') * pow(jump(phi_n),2) * dS) + 0.25 * pow(pow(phi_n,2)-1,2)*dx)
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

pic = plot(phi)
plt.title("Función de campo de fase en t = %.4f" %(t))
plt.colorbar(pic)
plt.show()


plt.plot(np.linspace(0,T,num_steps+1),E, color='red')
plt.title("Energía discreta")
plt.xlabel("Tiempo")
plt.ylabel("Energía")
if(savepic):
    plt.savefig("fig/DG-Eyre_nt-%d_energia.png" %(num_steps))
plt.show()


# "Main function"

if ( __name__ == '__main__' ):

    # Create mesh and define function space
    nx = ny = 100 # Boundary points
    print("nx = ny = %d" %(nx))
    mesh = UnitSquareMesh(nx,ny)

    # Random initial data
    degree = 1
    random.seed(1)
    class Init_u(UserExpression):
        def eval(self, values, x):
            values[0] = random.uniform(-0.49,0.51)

    phi_0 = Init_u(degree=degree) # Random values between -0.01 and 0.01

    T = 0.05            # final time
    num_steps = 100     # number of time steps
    dt = T / num_steps # time step size
    t = 0

    print("Iteraciones:")

    for i in range(num_steps):

        print("\nIteración %d:" %(i+1))

        # Update current time
        t += dt
        CH_time_step(mesh, degree, t, dt, phi_0)
