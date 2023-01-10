from dolfin import *
import matplotlib.pylab as plt

parameters['linear_algebra_backend'] = 'PETSc'
parameters["mesh_partitioner"] = "ParMETIS" # "SCOTCH"
# parameters["num_threads"] = 1;

p = 1
n = 16
mesh = UnitSquareMesh(n, n)


def refine_cell(cell):
    p = cell.midpoint()
    if p.distance(origin) < 0.1:
        do_refine = True
    else:
        do_refine = False
    return do_refine


cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
cell_markers.set_all(False)
origin = Point([0.0, 0.0])
for cell in cells(mesh):
    cell_markers[cell] = refine_cell(cell)


mesh = refine(mesh, cell_markers)
plot(mesh)
plt.show()

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
