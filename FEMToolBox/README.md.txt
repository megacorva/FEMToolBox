# FEMToolBox

FEMToolBox is an educational MATLAB toolbox for studying the Finite Element Method (FEM) on two-dimensional polygonal domains.

It provides accessible data structures and numerical routines that let users choose the computational mesh, inspect the finite element representation, and control the numerical scheme used to solve partial differential equations.

## Motivation

MATLAB’s built-in PDE functions can generate a mesh, construct a numerical scheme, and solve the resulting system through a largely black-box workflow. This is convenient for applications, but it can hide important details from users studying the Finite Element Method.

FEMToolBox exposes these components so that users can:

- provide or construct their own triangular mesh;
- inspect nodes, elements, boundary nodes, and internal nodes;
- assemble and examine finite element matrices;
- refine meshes and finite element functions;
- select coefficients and time-stepping parameters;
- study how the numerical method produces an approximate solution;
- visualize spatial and time-dependent finite element solutions.

The project is intended primarily for academic study, experimentation, and teaching rather than as a replacement for MATLAB’s production PDE solvers.

## Scope

The toolbox currently focuses on:

- two-dimensional polygonal domains;
- conforming triangular meshes;
- continuous piecewise-linear (`P1`) finite elements;
- homogeneous Dirichlet boundary conditions;
- elliptic and parabolic partial differential equations;
- obstacle problems;
- uniform spatial refinement;
- interpolation and truncation-error estimates.

## Main Components

### `P1Mesh`

Represents a two-dimensional triangular mesh.

A mesh stores:

- node coordinates;
- triangle connectivity;
- boundary-node indices;
- internal-node indices;
- the detected convexity of the domain.

It also supports:

- mesh visualization;
- uniform refinement, dividing every triangle into four triangles;
- assembly of diffusion, convection, and reaction matrices;
- computation of interpolation-error constants;
- computation of the domain diameter.

### `FEFunc`

Represents a spatial finite element function through its nodal values on a `P1Mesh`.

It supports:

- visualization of a finite element solution;
- uniform refinement while preserving the represented function.

### `TFEFunc`

Represents a time-dependent finite element function.

Its nodal values are stored at a sequence of sampling times. It supports:

- animated visualization;
- exporting an animation as an MPEG-4 video;
- uniform spatial refinement;
- uniform temporal refinement.