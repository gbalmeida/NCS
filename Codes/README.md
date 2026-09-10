# Source Code - Algorithms and Procedures

This directory contains the C++ implementation of the algorithms and cutting plane procedures proposed in the manuscript.

---

## 1. Requirements

- **Operating System:** Linux (tested on Ubuntu 24.04 LTS, 64-bit)
- **Compiler:** `g++` (GCC) with 64-bit support (`-m64`)
- **Mathematical Programming Solver:** IBM ILOG CPLEX Optimization Studio (version 22.1.1 or compatible)

---

## 2. Configuring the Makefile

Before compiling, open the provided `Makefile` in any text editor and adjust the `ILOG` variable (line 3) to point to your local CPLEX Optimization Studio installation directory:

```makefile
ILOG=/opt/ibm/ILOG/CPLEX_Studio2211
```

---

## 3. Compilation

To build the executable, open a terminal inside this directory and run:

```bash
make
```

to generate the binary executable named `cs`.

---

## 4. Execution

Once compiled, run the algorithm directly from the terminal by passing the path to the instance file as an argument:

```bash
./cs <path_to_instance>/instance
```

### Example:
Execution example for instance `i300_1.plc` when the file is located in the same directory as the executable `cs`:

```bash
./cs i300_1.plc
```


