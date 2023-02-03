**Dependencies:** numpy, mergedeep, mpi4py, matplotlib, pyyaml
### Exit codes:
**Exit code 0:** Process terminated normally.  
**Exit code 11:** Process terminated due to either a shuffle configuration error or NaN was encountered in the likelihoods. The process will be restarted automatically.  
**Exit code 12:** Process terminated as NestOR ran out of maximum allowed iterations. The process will not be restarted.  
**Exit code 13:** Process terminated due to *Math domain error* in analytical uncertainty calculation. This happened probably because the run terminated too early.  