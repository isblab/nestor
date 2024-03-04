Place the private header files in this directory. They will be
available to your code with

     #include <IMP/nestor/internal/myheader.h>

All headers should include `IMP/nestor/nestor_config.h` as their
first include and surround all code with `IMPNESTOR_BEGIN_INTERNAL_NAMESPACE`
and `IMPNESTOR_END_INTERNAL_NAMESPACE` to put it in the
IMP::nestor::internal namespace and manage compiler warnings.
