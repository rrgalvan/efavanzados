from dolfin import *

parameters['linear_algebra_backend'] = 'PETSc'
parameters["mesh_partitioner"] = "ParMETIS" # "SCOTCH"
# parameters["num_threads"] = 1;

p = 1
n = 1024 #640 #320 #64
mesh = UnitSquareMesh(n, n)
Q = FunctionSpace(mesh, "CG", p)
v = TestFunction(Q)
u = TrialFunction(Q)
a = dot(grad(u), grad(v))*dx
L = Constant(1.0)*v*dx

def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS
bc = DirichletBC(Q, Constant(0.0), boundary)

w = Function(Q)
solve(a==L, w, bc)
print("Vale")
