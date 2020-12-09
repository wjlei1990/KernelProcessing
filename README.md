# KernelProcessing
Kernel Processing program for Adjoint Tomography

0. Remove source mask (optional)

1. Sum event kernels to one kernel file

2. smooth each kernel to single files

3. merge smoothed kernels into one single file

4. compute hessian kernels for vp, vs and eta (optional)

5. precondition kernels

6. calculate search direction
  0. Steepest Descent
  1. Conjuate Gradient
  2. L-BFGS

7. update models using the search direction
