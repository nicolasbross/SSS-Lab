## Overview
The goal of this framework provides is to provide a rigorous numerical environment for analyzing elliptic operators of the form:
$$-\Delta u + V(x)u = \lambda u$$

## Key Features

### 1. 
* **Operator Support:** Handles a diverse class of potentials $V(x)$:
    * Analytic: Free, Harmonic (Confining), Coulomb, Radial.
    * Singular/Rough: Magnetic potentials, and potentials with local singularities.
* **Quantitative information:** 
* $\|u\|_{W^{k,p}}$ and $[u]_{C^{0,\alpha}}$ estimates.
* Sampling of frequency functions and doubling indices for unique continuation quantification.
* **Localization Measures:** Entropy-based localization and spectral gaps.

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



