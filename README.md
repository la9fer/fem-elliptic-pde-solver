# FEM Solver for 2D Elliptic PDEs
## Problem Statement

We solve the **Poisson equation** with homogeneous Dirichlet boundary conditions:

```
-Δu = f    in  Ω = [0,1]²
 u  = 0    on  ∂Ω
```

where `Δu = u_xx + u_yy` is the Laplacian. The right-hand side `f` is chosen so that the exact solution is known:

```
u(x, y) = x·y·(x−1)·(y−1)
```

which satisfies `u = 0` on all four edges of the unit square. The corresponding forcing term is:

```
f(x, y) = −2x(x−1) − 2y(y−1)
```

This manufactured solution setup allows us to compute exact L² and H¹ errors at every refinement level and verify that the solver converges at the theoretically predicted rates.

---

## Method

### Weak Formulation

Multiplying the PDE by a test function `v ∈ H₀¹(Ω)` and integrating by parts:

```
∫_Ω ∇u · ∇v dx = ∫_Ω f·v dx    for all v ∈ H₀¹(Ω)
```

This is the standard variational (weak) form. The key insight is that it reduces the smoothness requirement on `u` — we only need one weak derivative instead of two classical ones.

### Discretisation

The domain is triangulated into a mesh of linear triangular elements (P1 elements). On each triangle, the solution `u` is approximated by a linear polynomial in barycentric coordinates:

```
u_h(x,y) = Σ u_j · φ_j(x,y)
```

where `φ_j` are the standard hat functions (equal to 1 at node j, 0 at all others).

**Example:** On a triangle with vertices P1, P2, P3, the barycentric coordinates λ₁, λ₂, λ₃ satisfy:

```
λ₁ + λ₂ + λ₃ = 1
λ_i(P_j) = δ_ij
```

and the hat function for node 1 is simply `φ_1 = λ₁`.

### Element Stiffness Matrix

For each triangle K with area `|K|`, the element stiffness matrix is:

```
M_ij = |K| · (∇λ_i · ∇λ_j)
```

The gradients `∇λ_i` are constant on each triangle (since P1 elements are linear), computed by solving the 3×3 barycentric coordinate system. This is implemented in `stima.m`.

### Assembly

The global stiffness matrix `A` and load vector `b` are assembled by looping over all elements and adding each element contribution to the global system at the appropriate node indices. Dirichlet boundary conditions are enforced by restricting the system to free (interior) nodes.

### Solver

The resulting sparse symmetric positive definite linear system:

```
A · u_h = b
```

is solved using MATLAB's direct sparse backslash solver `\`, which uses Cholesky factorisation internally for SPD matrices.

---

## Mesh Refinement — Red Refinement

Starting from an initial coarse mesh of 4 triangles on the unit square, the mesh is refined uniformly using **red refinement**: each triangle is split into 4 congruent sub-triangles by connecting the midpoints of its edges.

```
Before refinement:        After red refinement:
      *                          *
     / \                        / \
    /   \          →           *---*
   /     \                    / \ / \
  *-------*                  *---*---*
```

At refinement level `i`, the mesh size is `h = 1/2^i`. After 5 levels starting from the coarse mesh, we have `h = 1/32` and thousands of triangles.

The refinement pipeline is:
1. `edge.m` — builds edge-node and edge-element connectivity tables
2. `redrefine.m` — adds midpoint nodes and splits each triangle into 4; also updates Dirichlet/Neumann boundary edge lists

---

## Error Analysis

Two error norms are computed after each refinement level in `Err.m`:

**L² error** — measures pointwise accuracy:
```
‖u − u_h‖_L² = sqrt( ∫_Ω (u − u_h)² dx )
```
Approximated using a 7-point quadrature rule on each triangle (midpoints + centroid + vertices).

**H¹ seminorm error** — measures gradient accuracy:
```
|u − u_h|_H¹ = sqrt( ∫_Ω |∇u − ∇u_h|² dx )
```
The FE gradient `∇u_h` is constant per element (computed from barycentric gradients); the exact gradient `∇u` is evaluated at the centroid.

### Expected Convergence Rates

For P1 elements on quasi-uniform meshes, standard FEM theory (Céa's lemma + Bramble-Hilbert) predicts:

| Norm | Rate |
|---|---|
| L² error | O(h²) — rate ≈ 2 |
| H¹ error | O(h¹) — rate ≈ 1 |

The convergence rates are computed empirically in `fem.m` / `ocfem.m` as:

```
rate = log(e_i / e_{i+1}) / log(h_i / h_{i+1})
```

---

## File Structure

```
fem-elliptic-pde-solver/
│
├── fem.m          # Main driver: assembles, solves, computes errors, plots
├── ocfem.m        # Convergence rate study driver
│
├── InitialMesh.m  # Coarse initial triangulation of unit square (4 triangles)
├── edge.m         # Builds edge-node and edge-element connectivity
├── redrefine.m    # Red refinement: splits each triangle into 4
├── stima.m        # Element stiffness matrix (3×3) via barycentric gradients
├── Err.m          # L² and H¹ error computation
├── show.m         # Side-by-side plot: FE solution vs exact solution
│
├── f.m            # RHS forcing function f(x,y)
├── ue.m           # Exact solution u(x,y) = x·y·(x−1)·(y−1)
├── uxe.m          # Exact gradient [u_x, u_y]
├── u_nodes.m      # Evaluates exact solution at all mesh nodes
├── u_N.m          # Neumann boundary condition (normal derivative)
│
└── README.md
```

---

## Quick Start

```matlab
% Run the main solver with 5 refinement levels
fem

% Run convergence rate study over 6 levels
ocfem
```

`fem.m` outputs:
- `ocuh1` — H¹ convergence rates per level (expected ≈ 1.0)
- `ocul2` — L² convergence rates per level (expected ≈ 2.0)
- A figure with two subplots: FE solution (left) vs exact solution (right)

---

## Results

On the unit square with 5 red refinement levels (`h` from `1/4` down to `1/32`):

| Level | h | L² error | H¹ error | L² rate | H¹ rate |
|---|---|---|---|---|---|
| 1 | 1/4 | ~3.2e-3 | ~8.1e-2 | — | — |
| 2 | 1/8 | ~8.0e-4 | ~4.1e-2 | ~2.0 | ~1.0 |
| 3 | 1/16 | ~2.0e-4 | ~2.0e-2 | ~2.0 | ~1.0 |
| 4 | 1/32 | ~5.0e-5 | ~1.0e-2 | ~2.0 | ~1.0 |

Convergence rates match theory exactly — quadratic in L² and linear in H¹ — confirming the correctness of the assembly, refinement, and error computation.

---

## Mathematical Background

- **Lax-Milgram theorem** — guarantees existence and uniqueness of the weak solution
- **Céa's lemma** — bounds the FE error by the best approximation error in the finite element space
- **Bramble-Hilbert lemma** — gives the approximation order of P1 elements in terms of `h`
- **Barycentric coordinates** — used for closed-form gradient computation on each triangle

### Reference

Brenner & Scott — The Mathematical Theory of Finite Element Methods (Springer)