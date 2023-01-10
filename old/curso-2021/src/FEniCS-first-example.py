from fenics import *
import matplotlib.pylab as plt

#
# 1) Mesh
#
n_intervals = 10
a, b = 0, 1
mesh = UnitSquareMesh(n_intervals, n_intervals)
plot(mesh)
plt.show()

#
# 2) Function spaces
#
polynom_degree = 2
U = FunctionSpace(mesh, "Lagrange", polynom_degree)
u = TrialFunction(U)
v = TestFunction(U)

#
# 3) Variational fomulation
#
def boundary(x):
    return ( x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS
             or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS)
bc0 = DirichletBC(U, 0, boundary)

f = Expression("sin(pi*x[0])*sin(pi*x[1])", degree=1)
a = dot(grad(u),grad(v)) * dx
b = f*v * dx

#
# 4) Problem solving and plotting solution
#
sol = Function(U)
solve(a == b, sol, bc0)

plot(sol)
plt.show()
