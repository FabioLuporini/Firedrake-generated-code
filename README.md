# Firedrake-generated-code
Finite element kernels generated through Firedrake

How to interpret the generated code
===================================
Each .c file in /toms includes the generated code for a specific test case.  In
a file there are always two functions: the actual assembly kernel (at the top of
the file) and the so-called "wrapper function", which invokes the kernel for
each mesh element. PyOP2, the software module used by Firedrake to apply kernels
to mesh elements, jumps from Python- to C-land by calling the wrapper function.

The assembly kernels in /toms have been generated through various compilers:
FFC ("original"), FFC's optimized quadrature mode ("quad"), FFC's tensor
contraction mode ("tensor"), UFLACS ("uflacs"), and COFFEE ("coffee"). A
filename has the following format: ::

    problem_mode_nfX_qY.c

Where ``problem`` usually identifies the operator being assembled (e.g.,
``mass`` for the mass matrix); ``mode`` is one of the aforementioned compilers
(e.g., ``original``, ``tensor``); ``nf`` indicates the number of pre-multiplying
functions (``X`` in total); ``q`` indicates the polynomial degree (of value
``Y``).
