Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/nestor/myheader.h>

All headers should include `IMP/nestor/nestor_config.h` as their
first include and surround all code with `IMPNESTOR_BEGIN_NAMESPACE`
and `IMPNESTOR_END_NAMESPACE` to put it in the IMP::nestor namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
