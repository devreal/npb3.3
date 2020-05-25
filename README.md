# npb3.3
Ports of the popular NAS Parallel Benchmark BT-MZ to different task-based programming models.

- `BT-MZ-CXX`: C++ port of the Fortran reference implementation using MPI and OpenMP worksharing constructs
- `BT-MZ-CXX-omptasks`: variant of the C++ port using tasks instead of parallel loops, fork-join model with MPI outside of tasks
- `BT-MZ-CXX-omptasks-detach-merged`: variant of the C++ port using OpenMP tasks, using detached tasks for MPI communication (with send/recv *merged* into one task). This port relies on MPI Continuations, an experimental extension to MPI for which a prototypical implementation is available at https://github.com/devreal/ompi/tree/mpi-continue-master.
- `BT-MZ-CXX-ompss2tasks-tampi-events-merged`: C++ port using OmpSs-2 tasks and TAMPI with its interface using OmpSs-2 *events* (with send/recv *merged* into one task)
- `BT-MZ-CXX-ompss2tasks-continue-events-merged`: C++ port using OmpSs-2 tasks with the *events* interface to block the release of dependencies while communication is outstanding (with send/recv *merged* into one task). This port relies on MPI Continuations, an experimental extension to MPI for which a prototypical implementation is available at https://github.com/devreal/ompi/tree/mpi-continue-master.

The C++ port uses a generic N-dimensional array abstraction (`NArray.h`) that requires C++17 features. Unfortunately, the OmpSs-2 compiler only supports up to C++11 so the computational code was split from the OmpSs statements. The OmpSs parts are found in the respective `*_omp.cc` files, the parts using the NArray class is found in `*_impl.cc` and a header glues them together.
