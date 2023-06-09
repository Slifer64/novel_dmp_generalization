# Experiments

The experiments were carried out on a desktop PC with an Intel® Core™ i7-9700 processor.
The code was written in c++ and organized in ros packages.
Since there are a lot 3rd party software package dependencies, and inter-dependencies between our own custom packages, we provide only part of the source code for the experiments and the DMP library. If needed, a docker image could probably be added in the near future.

---

## Code strucute

- The DMP++ described in the paper is implemented in `src/DMP_pp.h` and `src/DMP_pp.cpp`. It uses the basic functionality of our DMP library implementation, organized as a ros package in `src/lib/gmp_lib`. The latter depends also on the (customized) OSQP library `src/lib/osqp_lib`, which is however not used for the implementation of DMP++.
- The DMP models used for pick and place are stored in `src/pick_model.bin` and `src/place_model.bin`. The control parameters used are included in `src/ctrl_params.yaml`.
- The controller is implemented in `src/controller.cpp`. 
The implementation of some auxiliary class objects and functions is not incorporated. One could adopt the same structure-logic with their own implementation of these auxiliary utilities.


