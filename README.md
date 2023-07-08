# PS-SBFEM
**Title:** A Novel Solution for Seepage Problems Implemented in the Abaqus UEL Based on the Polygonal Scaled Boundary Finite Element Method

**Abstract:** 
The scaled boundary finite element method (SBFEM) is a semianalytical computational scheme based on the characteristics of the
finite element method (FEM) and boundary element method that combines their respective advantages. In this paper, the SBFEM
and polygonal mesh technique are integrated into a new approach to solve steady-state and transient seepage problems. The
proposed method is implemented in Abaqus employing a user-defined element (UEL). A detailed implementation of the
procedure is presented in which the UEL element is defined, the internal variables RHS and AMATRX are updated, and the
stiffness/mass matrix is solved using eigenvalue decomposition. Several benchmark problems are solved to verify the proposed
implementation. The results show that the polygonal element of the polygonal SBFEM (PSBFEM) is more accurate than the
standard FEM element of the same element size. For transient problems, the results for the PSBFEM and FEM are in excellent
agreement. Hence, the proposed method is robust and accurate for solving steady-state and transient seepage problems.

# How to link the MKL
- We developed a UEL in Abaqus 6.14. In order to solve for eigenvalues and eigenvectors, we used the MKL library. Therefore, before using this UEL, it is necessary to link the MKL library.
- Fortunately, Abaqus comes with the MKL library, so we just need to modify the environment file to enable the linking. Please refer to the "abaqus_v6.env" file for specific instructions on how to do this. 

# How to use the UEL
- You can run UEL with the command
```dos
call abaqus job=xxx user=sbfem2DSeepage.for int 
```

# Citation
If you use PSBFEM-Abaqus for academic research, you are encouraged to cite the following paper:

```
@article{yang2022novel,
  title={A Novel Solution for Seepage Problems Implemented in the Abaqus UEL Based on the Polygonal Scaled Boundary Finite Element Method},
  author={Yang, Yang and Zhang, Zongliang and Feng, Yelin and Wang, Kun and others},
  journal={Geofluids},
  volume={2022},
  year={2022},
  publisher={Hindawi}
}
```
