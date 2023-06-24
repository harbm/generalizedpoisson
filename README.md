# GeneralizedPoissonSolver
This repo contains an efficient solver for the 3-D generalized Poisson Equation. 

# Installation
- Ensure you have MATLAB installed on your workstation. It should work with most versions but I've only tested it with 2022a, 2022b, and 2023a. 
- Clone or download the repo.
- 

## FAQ

### What's the difference between the Poisson Equation and the _generalized_ Poisson Equation?   


<b>Short Answer</b>: The Poisson Equation is applicable to systems where the dielectric permitivity is _constant_. The _generalized_ Poisson equation is further applicable to systems where the dielectric permitivity is _spatially varying_.

<b>Long Answer</b>: Let's start from the most general formulation of Gauss's law.

$$
\nabla \cdot \mathbf{D}({\mathbf{r}}) = \rho({\mathbf{r}})
$$

Here, ${\bf{r}}= x{\bf{\hat{x}}} + y{\bf{\hat{y}}} +z{\bf{\hat{z}}}$, $ \mathbf{D}$ is the electric displacement field, and $\rho$ is the electric charge density. 

In a _linear_ material (i.e. one where the polarization depends _linearly_ on the electric field) we can use the relation

$$
\mathbf{D}({\mathbf{r}}) = \varepsilon_0 \varepsilon_r({\mathbf{r}}) \  \mathbf{E}({\mathbf{r}})
$$

to rewrite Gauss's law in terms of the electric field $\mathbf{E}$

$$
\nabla \cdot \left[\varepsilon_0\varepsilon_r({\mathbf{r}}) \mathbf{E}({\mathbf{r}})\right] = \rho({\mathbf{r}})
$$}

where $\varepsilon_0$ is the permitivity of free-space and $\varepsilon_r$ is the relative permativity. Further rewriting Gauss's law in terms of the electric potential function $\phi$ by making the substitution $E({\mathbf{r}})=-\nabla \phi({\mathbf{r}})$, we arrive at an expression for the _generalized Poisson Equation_:

$$
\nabla \cdot \left[\varepsilon_r({\mathbf{r}}) \nabla \phi({\mathbf{r}})\right] = -\frac{\rho({\mathbf{r})}}{\varepsilon_0}
$$

To arrive at the more familiar expression for the _Poisson Equation_, we further assume a _uniform_ dielectric permitivity $\varepsilon_r({\mathbf{r}})=\varepsilon_r$ 

$$
\nabla^2\phi({\mathbf{r}}) = - \frac{\rho({\mathbf{r})}}{\varepsilon_0\varepsilon_r}
$$
## What are the boundary conditions supported by this solver?
Dirichlet and Neumann boundary conditions are both supported, as well as combinations of both (e.g. Dirichlet on some sides and Neumann on other sides.)

### Different PDEs

### What license is this software distributed under?
This software is distributed under the terms of the MIT License. 

### How should I cite this software? 

If you utilize this software in academic research, I request that you cite it appropriately. Please include the following in your citation:

    Mohammed Harb. (2022). GeneralizedPoissonSolver. 1.0. Retrieved from https://github.com/harbm/generalizedpoisson

Or in BibTeX format: 



    @misc{AuthorYear,
        author = {Your Name},
        title = {SoftwareName},
        year = {Year},
        publisher = {GitHub},
        journal = {GitHub repository},
        howpublished = {\url{https://github.com/username/repository}}
    }