# Optimal Shape Design for Poisson Equation

We study the shape design problem through the minimization of the cost function

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BJ%7D%20%5Cleft%28%20%5Ctheta%2C%20y%20%5Cright%29%20%3D%20%5Cfrac%7B1%7D%7B2%7D%20%5Cint_%7B%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%7D%20%5Cleft%28%20y%20-%20y_d%20%5Cright%29%20%5E2%20%5C%2C%20%5Cmathrm%7Bd%7D%20%5COmega%2C">
</p>

where <img src="https://latex.codecogs.com/gif.latex?y%20%5Cin%20L%5E2%20%5Cleft%28%20%5COmega%20%5Cright%29"> is the state variable, <img src="https://latex.codecogs.com/gif.latex?y_d%20%5Cin%20L%5E2%20%5Cleft%28%20%5COmega%20%5Cright%29"> is a target function, and <img src="https://latex.codecogs.com/gif.latex?%5Ctheta%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29"> the normal displacement to a reference boundary,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma%20%5Cleft%28%20%5Ctheta%20%5Cright%29%20%3D%20%5Cleft%5C%7B%20%5Cmathbf%7Bx%7D%20&plus;%20%5Ctheta%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cmathbf%7Bn%7D%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%3A%20%5Cmathbf%7Bx%7D%20%5Cin%20%5CGamma_%7B0%7D%20%5Cright%5C%7D.">
</p>

The problem is subject to the following elliptic PDE in <img src="https://latex.codecogs.com/gif.latex?%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%5Csubset%20%5Cmathbb%7BR%7D%5Ed"> with Dirichlet boundary conditions on <img src="https://latex.codecogs.com/gif.latex?%5CGamma%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%3D%20%5CGamma_w%20%5Ccup%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29">,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20-%5CDelta%20y%20%3D%20f%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%5C%5C%20y%20%3D%20y_%7Bw%7D%20%26%20%5Ctext%7Bin%20%7D%20%5CGamma_w%2C%5C%5C%20y%20%3D%20y_%7Bs%7D%20%26%20%5Ctext%7Bin%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5Cend%7Bcases%7D">
</p>

with <img src="https://latex.codecogs.com/gif.latex?f%20%5Cin%20L%5E2%20%5Cleft%28%20%5COmega%20%5Cright%29">.

We consider the set of admissible domains whose measure is fixed,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BU%7D_%7Bad%7D%20%3D%20%5Cleft%5C%7B%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%3A%20%5Clvert%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%5Crvert%20%3D%20%5COmega_0%20%5Cright%5C%7D%2C">
</p>

and aim at finding the optimal shape that minimizes the cost function. In order to do so, we use the steepest descent method with descent direction given by

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cdelta%20%5Ctheta%20%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%3D%20-%20%5Cleft%5B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cleft%28%20y%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20-%20y_d%20%5Cright%29%20%5E2%20&plus;%20%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cfrac%7B%5Cpartial%20y%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20&plus;%20q%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cright%5D%2C">
</p>

with <img src="https://latex.codecogs.com/gif.latex?v"> solution of the adjoint problem

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20-%20%5CDelta%20v%20%3D%20y%20-%20y_d%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20v%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%3D%20%5CGamma_w%20%5Ccup%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

The normal displacement is updated at every iteration,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Ctheta%5E%7B%5Cleft%28%20n%20&plus;%201%20%5Cright%29%7D%20%3D%20%5Ctheta%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20&plus;%20%5Cepsilon%20%5Cdelta%7B%5Ctheta%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%7D%2C">
</p>

with <img src="https://latex.codecogs.com/gif.latex?\epsilon"> sufficiently small. The Lagrange multiplier is computed in order to ensure that the volume contraint is fulfilled, thus

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?q%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%3D%20-%20%5Cfrac%7B1%7D%7B%5Cleft%7C%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%5Cright%7C%7D%20%5Cint_%7B%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%7D%20%5Cleft%5B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cleft%28%20y%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20-%20y_d%20%5Cright%29%20%5E2%20&plus;%20%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cfrac%7B%5Cpartial%20y%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cright%5D%20%5Cmathrm%7Bd%7D%20%5CGamma.">
</p>

## Mesh Motion Solver

The solutions to the primal and adjoint problems are commonly approximated by means of numerical methods, such as the finite element method (FEM) or the finite volume method (FVM). In order to apply these techniques, the domain under study must be tessellated with a mesh. Nevertheless, applying the displacement directly to the controlled boundary will deteriorate the surrounding elements after a few iterations and the computation will crash if the interior nodes of the domain are not reallocated. In order to avoid this, the domain can be re-meshed after a number of iterations. However, this can be very expensive as a completely new mesh must be generated. A commonly used alternative is to move the interior nodes of the mesh according to the displacements prescribed on the boundary.  By doing this the number of elements and the nodes connectivities remain the same, only the mesh nodes positions are updated. The Solid Body Rotation Stress method and the Laplacian smoothing have been used in OpenFOAM.

### Linear Elasticity

Mesh motion can be achieved by treating the mesh as an elastic body and solving the equations of Linear Elasticity for solids with prescribed displacements on the domain boundary,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20%5Cpartial_j%20%5Cleft%28%20G%20%5Cpartial_j%20u_i%20%5Cright%29%20&plus;%20%5Cpartial_j%20%5Cleft%28%20G%20%5Cpartial_i%20u_j%20%5Cright%29%20&plus;%20%5Cpartial_i%20%5Cleft%28%20%5Clambda%20%5Cpartial_j%20u_j%20%5Cright%29%20%3D%200%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20u_i%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_w%2C%20%5C%5C%20u_i%20%3D%20%5Cdelta%20%5Ctheta%20n_i%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

### Solid Body Rotation (SBR) Stress

The Linear Elasticity model fails when dealing with rotating meshes. This can be mitigated by selecting the material properties in a manner,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%202%20%5Cpartial_j%20%5Cleft%28%20%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cpartial_j%20u_i%20%5Cright%29%20&plus;%20%5Cpartial_j%20%5Cleft%5B%20%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cleft%28%20%5Cpartial_i%20u_j%20-%20%5Cpartial_j%20u_i%20-%20%5Cdelta_%7Bij%7D%20%5Cpartial_k%20u_k%20%5Cright%29%20%5Cright%5D%20%3D%200%2C%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20u_i%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_w%2C%20%5C%5C%20u_i%20%3D%20%5Cdelta%20%5Ctheta%20n_i%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

### Laplacian Smoothing

A simple and widely used practice is to solve a Laplace equation with the prescribed displacements as boundary conditions,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20%5CDelta%20u_i%20%3D%200%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20u_i%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_w%2C%20%5C%5C%20u_i%20%3D%20%5Cdelta%20%5Ctheta%20n_i%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

The behavior of the Laplacian smoothing can be improved by adding a non-uniform diffusivity term,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cpartial_k%20%5Cleft%28%20%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cpartial_k%20u_i%20%5Cright%29%20%3D%200%20%5Cquad%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.">
</p>

The diffusion field <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%3A%20%5Cmathbb%7BR%7D%5Ed%20%5Crightarrow%20%5Cmathbb%7BR%7D"> decreases with the distance to the controlled boundary. A common choice is to make it depend on the inverse of the distance as

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%3D%20%5Cfrac%7B1%7D%7Bd%5Em%7D%2C">
</p>

so that the nodes next to the deforming boundaries move with similar displacements as those on the boundary.

## Getting Started

The solver must be compiled in the terminal. It is advisable to first clean previous compilations with

```
wclean
```

and then use

```
wmake
```

### Prerequisites

OpenFOAM C++ library must be installed in order to compile the code.

The OpenFOAM distribution provided by the [OpenFOAM Foundation](https://openfoam.org/) was used.

## Running a Case

In order to run the solver move to the case folder _poissonCGAdjoinFoamCase_ and type in the command line

```
./Allprepare

poissonCGAdjointFoam
```

The shape optimization method described in the previous section has been tested with a simple example. The Poisson equation is posed in a two-dimensional circular domain with boundary

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma_w%20%3D%20%5Cleft%5C%7B%20%5Cleft%28%20x%2C%20y%20%5Cright%29%20%5Cin%20%5Cmathbb%7BR%7D%5E2%3A%20%5Csqrt%7Bx%5E2%20&plus;%20y%5E2%7D%20%3D%20R_w%20%5Cright%5C%7D.">
</p>

The reference geometry to be optimized is an inner hole with boundary given by

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma_s%20%5Cleft%28%200%5Cright%29%20%3D%20%5Cleft%5C%7B%20%5Cleft%28%20x%2C%20y%20%5Cright%29%20%5Cin%20%5Cmathbb%7BR%7D%5E2%20%3A%20%5Csqrt%7B%5Cleft%28%20x%20-%20c_x%20%5Cright%29%5E2%20&plus;%20%5Cleft%28%20y%20-%20c_y%20%5Cright%29%5E2%7D%20%3D%20R_s%20%5Cright%5C%7D.">
</p>

<p align="center">
  <img src="poissonOptShapeFoamCase/figs/fig1.png" width="400" height="300">
</p>


The _poissonCGAdjointFoam_ solver has been tested in a square domain <img src="https://latex.codecogs.com/gif.latex?%5B0%2C%201%5D%20%5Ctimes%20%5B0%2C%201%5D"> with zero Dirichlet boundary conditions and <img src="https://latex.codecogs.com/gif.latex?%5Cbeta%20%3D%2010%5E%7B-3%7D%2C10%5E%7B-4%7D%2C10%5E%7B-5%7D%2C10%5E%7B-6%7D">. The target function is <img src="https://latex.codecogs.com/gif.latex?y_d%20%3D%20xy%20%5Csin%20%5Cleft%28%20%5Cpi%20x%20%5Cright%29%20%5Csin%20%5Cleft%28%20%5Cpi%20y%20%5Cright%29">.

<p align="center">
  <img src="poissonCGAdjointFoamCase/cg_J.png">
</p>

<p align="center">
  <img src="poissonCGAdjointFoamCase/cg_Jy.png">
</p>

### Warning

It might be needed to use 

```
sed -i -e 's/\r$//' filename
```

and

```
chmod +x filename
```

in order to be able to execute 

```
./filename
```

## Author

* **Jose Lorenzo Gomez**

## Acknowledgments

This project has received funding from the European Research Council (ERC) under the European  Union’s Horizon 2020 research and innovation programme (grant agreement No. 694126-DyCon).
 
[DyCon Webpage](http://cmc.deusto.eus/dycon/)

## References
