# P1 Finite Element Solver for the 2D Poisson Equation

A from-scratch implementation of the **linear finite element method (P1 FEM)** for solving the two-dimensional Poisson equation on the unit square. The project includes triangular mesh generation, element stiffness matrix construction, sparse global assembly, uniform red refinement, and numerical error estimation.

Implemented in MATLAB

## Problem formulation

The solver considers the Poisson boundary-value problem

\[
-\Delta u=f \quad \text{in } \Omega=(0,1)^2,
\]

\[
u=0 \quad \text{on } \partial\Omega.
\]

A manufactured exact solution is used to evaluate the numerical error:

\[
u(x,y)=xy(x-1)(y-1).
\]

The corresponding forcing function is

\[
f(x,y)=-2x(x-1)-2y(y-1).
\]

The exact solution satisfies the homogeneous Dirichlet boundary condition.

## Numerical method

The problem is discretised using continuous, piecewise-linear triangular finite elements.

### 1. Initial mesh and refinement

The initial mesh consists of four triangular elements covering the unit square, with an interior node at its centre.

Uniform red refinement divides each triangle into four smaller triangles by introducing edge midpoints.

The refinement procedure is implemented through:

- `InitialMesh.m`: constructs the initial triangulation and boundary edges.
- `edge.m`: identifies edge-node and edge-element connectivity.
- `redrefine.m`: creates new nodes, elements, and refined boundary edges.

### 2. Element stiffness matrix

For each triangular element, the local stiffness matrix is computed using the gradients of the linear barycentric basis functions.

The local matrix is assembled into the global sparse stiffness matrix.

Implementation: `stima.m`.

### 3. Global assembly and solution

The global stiffness matrix and load vector are assembled by iterating over the triangular elements.

The load vector uses a one-point quadrature rule evaluated at each element centroid.

The Dirichlet boundary values are imposed using the exact solution. The resulting sparse linear system is solved over the free degrees of freedom using MATLAB's sparse backslash operator.

Main implementations: `fem.m` and `ocfem.m`.

### 4. Error estimation

The numerical solution is compared with the manufactured exact solution.

Two error measures are computed:

- **L² error:** approximated using a seven-point quadrature rule per triangle, involving the vertices, edge midpoints, and centroid.
- **H¹-seminorm error:** compares the constant gradient of the P1 finite element solution with the exact gradient evaluated at the triangle centroid.

Implementation: `Err.m`.

## Expected convergence behaviour

For sufficiently regular solutions and quasi-uniform triangular meshes, standard P1 finite element theory predicts

\[
\|u-u_h\|_{L^2(\Omega)}=O(h^2),
\]

and

\[
|u-u_h|_{H^1(\Omega)}=O(h).
\]

The convergence study computes errors across successive uniform refinement levels and estimates the observed rates using

\[
r=\frac{\log(e_i/e_{i+1})}
        {\log(h_i/h_{i+1})}.
\]

Here, \(e_i\) denotes the error at refinement level \(i\), and \(h_i\) is the mesh-size indicator calculated in the convergence driver.

The expected asymptotic rates are approximately 2 in L² and 1 in the H¹ seminorm.

*Numerical convergence values should be added after reproducing the current source code results in MATLAB*

## Repository structure

| File | Description |
|---|---|
| `fem.m` | Main FEM driver: mesh refinement, assembly, solution, error computation, and visualisation |
| `ocfem.m` | Driver for the multi-level convergence study |
| `convergence_study.m` | Runs the convergence study and generates its summary plot |
| `InitialMesh.m` | Initial triangular mesh and boundary-edge construction |
| `edge.m` | Edge-node and edge-element connectivity |
| `redrefine.m` | Uniform red refinement of triangular elements |
| `stima.m` | Local element stiffness matrix |
| `Err.m` | L² and H¹-seminorm error estimation |
| `f.m` | Forcing function |
| `ue.m` | Manufactured exact solution |
| `uxe.m` | Exact gradient of the solution |
| `u_nodes.m` | Evaluation of the exact solution at mesh nodes |
| `u_N.m` | Neumann boundary data helper |
| `show.m` | Visual comparison of the finite element and exact solutions |

## Running the project

Place all MATLAB files in the same directory and run the following commands:

```matlab
% Run the FEM solver and visualise the solution
fem

% Run the multi-level convergence study
ocfem

% Generate the convergence summary and plot
convergence_study
```

The main solver and convergence driver use different numbers of refinement levels.

## Current scope and limitations

- The implementation uses P1 (piecewise-linear) triangular elements.
- The computational domain is the unit square.
- Mesh refinement is uniform rather than adaptive.
- The load vector uses centroid-based one-point quadrature.
- The H¹-seminorm error uses a centroid evaluation of the exact gradient.
- Global sparse matrix assembly is performed inside element loops, which may become computationally expensive on finer meshes.
- The implementation is intended as a numerical methods and FEM demonstration rather than a general-purpose PDE solver.

## References

- S. C. Brenner and L. R. Scott, *The Mathematical Theory of Finite Element Methods*, Springer
