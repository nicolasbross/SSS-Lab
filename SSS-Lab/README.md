## Overview
The goal of this framework is to provide a Terminal User Interface for analyzing elliptic operators of the form
$$-\Delta u + V(x)u = \lambda u$$. This includes treatment of the harmonic case $$-\Delta u = 0 $$ and the respective eigenvalue problem 
for the Laplacian. 

### 1. Features
* **Tools for Harmonic Analyiss** (constantly adding new):
   * Symbolic calculus for balls in $\mathbb{R}^n$, for integration over balls and spheres.
   * Computations of bases for spaces of spherical harmonics $\mathbb{R}^n$
   * Solutions of the Dirilecht problem over balls, ellipsoids and annular domains.
   * Computations of harmonic conjugates in $\mathbb{R}^2$.

* **Operator Support:** Handling of a diverse class of potentials $V(x)$:
    * Analytic: Free, Harmonic (Confining), Coulomb, Radial.
    * Singular/Rough: Magnetic potentials, and potentials with local singularities.
* **Quantitative information:** 
* Sampling of frequency functions and doubling indices.
* **Localization Measures:**
* Entropy-based localization and spectral gaps.

### 2. Design Goals
* Support explicit finite difference schemes and FEM-like weak formulations. Decouple the discretization layer from the solver logic, allowing for custom stencils.
* Implementation in modern **C++ (C++17/20)** utilizing the **Eigen** library for linear algebra backend.
* Develop an embedded **Lua** binding that allows for real-time parameter tuning (e.g., varying potential strength $\alpha$ in $V_\alpha$) and workflow orchestration without recompiling the C++ core.

## Project Structure

├── include/           # Header-only interfaces and template definitions
│   ├── operators/     # Abstract base classes for Elliptic Operators
│   ├── geometry/      # Grid and domain discretizations
│   └── analysis/      # Other related computations.
├── src/               # Compiled source implementations
│   ├── core/          # Solver logic and Eigen wrappers
│   └── binding/       # Lua bindings and API exposure
├── scripts/           # Lua and Bash scripts for experiment workflows
└── tests/             # Unit tests for convergence 



