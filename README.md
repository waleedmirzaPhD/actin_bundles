**Note:** This code reproduces **Figure 2b(I)** from the following scientific publication:  
> **"Theory of active self-organization of dense nematic structures in the actin cytoskeleton"**  
> *Waleed Mirza, Marco De Corato, Marco Pensalfini, Guillermo Vilanova, Alejandro Torres-Sànchez, Marino Arroyo*  
> [eLife Reviewed Preprint (2024)](https://elifesciences.org/reviewed-preprints/93097)




# Installation Instructions for HiperLife and Project Compilation

## 1. Install Required Dependencies
Ensure the following software and libraries are installed on your system:

- **C++ Compiler** (GCC, Clang, or MSVC)
- **MPI Library** (e.g., OpenMPI or MPICH)
- **Trilinos Library** (with Teuchos and Amesos enabled)
- **CMake** (for building the project)
- **Make** (or Ninja for faster compilation)

## 2. Install HiperLife

HiperLife is required for this project. Clone and install it from GitLab: https://gitlab.com/hiperlife/hiperlife/-/tree/release-4.3.0-branch


## 3. Clone and Compile HL Base Project
This project follows the HL Base Project structure. Follow the instructions from:

[HL Base Project](https://gitlab.com/hiperlife/hl-base-project)

To compile:

```sh
git clone https://gitlab.com/hiperlife/hl-base-project.git
cd hl-base-project
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

## 4. Build and Run the Application
Navigate to your project directory and follow these steps:

```sh
mkdir -p build && cd build
cmake ..
make -j$(nproc)
mpirun -np 4 ./your_executable input_file.config
```

## 5. Debugging and Verification
To check if MPI is correctly installed, run:

```sh
mpirun --version
```

To debug missing dependencies:

```sh
ldd your_executable | grep "not found"
```

For runtime debugging:

```sh
gdb --args mpirun -np 4 ./your_executable input_file.config
```

## 6. Cleaning Up
To remove compiled binaries and regenerate from scratch:

```sh
rm -rf build/
```

This document provides all necessary steps to install dependencies, compile, and run the project successfully.

