# Experiments

The experiments were carried out on a desktop PC with an Intel® Core™ i7-9700 processor.
The code was written in c++ and organized in ros packages.
Since there are a lot 3rd party software package dependencies, and inter-dependencies between our own custom packages, we provide only the source code for the experiments and the DMP library.

---

## DMP++

The DMP++ described in the paper is implemented in `src/DMP_pp.h` and `src/DMP_pp.cpp`, with basic documentation included in those files.
It uses the basic functionality of our DMP library implementation, organized as a ros package in `src/lib/gmp_lib`. The latter depends also on the (customized) OSQP library `src/lib/osqp_lib`, which is however not used for the implementation of DMP++.


---

## Execution
The core execution workflow is implemented in `src/controller.hpp`. Specifically:
- `executePick(bool reverse=false)`: For picking the carton from the conveyor belt.
- `executePlace(bool reverse=false)`: For placing the carton in the box.

`reverse=true` can be used to retract along the reverse path.
The purpose of these functions is to show the basic logic and control flow. Comments are included to facilitate understanding.
The implementation of some auxiliary class objects and functions is not incorporated. One could adopt the same structure-logic with their own implementation of these auxiliary utilities.