# MT-SBD-STM: Multi-Type Sparse blind deconvolution using the Riemannian Trust-Region Method (RTRM) on STM images

This is a MATLAB package that aims to deconvolute multi-type kernels and their corresponding activations

As sparse blind deconvolution is a nonconvex problem, using RTRM ensures that local minima will be found in the associated optimization objective.

This work is inspired by work in Single SBD-STM shown by [Cheung et al (2020)](https://www.nature.com/articles/s41467-020-14633-1), in which they formulated the Single SBD-STM problem: 
<img width="646" alt="image" src="https://github.com/user-attachments/assets/63946883-6cfa-44b0-b877-e28cfaf07c72" />

We extend the algorithm into multi type defects, as in most materials, we observe more than one type of defects, here is an illustration of the work: 
![image](https://github.com/user-attachments/assets/e48e51e4-87cc-4fa6-84c2-d0be148adbe3)

![image](https://github.com/user-attachments/assets/ddb843be-a8e6-4a6b-8871-60a3d1dd65ce)


