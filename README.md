[![PubMed](https://salilab.org/imp-systems/static/images/pubmed.png)](https://pubmed.ncbi.nlm.nih.gov/xyz/)     #TODO: Replace xyz with PubMed ID


# **NestOR: Nested Sampling-based Optimization of Representation for Integrative Structural Modeling**

## **Installation:**
### **Dependencies:**  
* IMP (compiled from the source code). See [IMP installation](https://github.com/salilab/imp)
* Python libraries: numpy, mergedeep, mpi4py, matplotlib, pyyaml, plotly
  
### **NestOR installation:**  
1. Compile IMP from the souce code to your choice of directory
2. Replace the `macros.py` in `imp/modules/pmi/pyext/src/` with the macros.py in the current repository. Make sure the file is named `macros.py` in the destination directory. 
3. Similarly, replace the restraints directory in `imp/modules/pmi/pyext/src/` with the restraints directory in the present repository.

## **Running NestOR:**

### Inputs 
1. Make the modeling script in the form as shown in the `example/modeling.py`. One will also need to make separate topology files for different candidate representations.
2. Set appropriate parameters in the `nestor_params.yaml` file. An example param file can be found in `examples/`.

### Run command
Run NestOR with the following command `python wrapper_vX.py $abs_nestor_params_path` where wrapper_vX.py is the wrapper script in this repository and `abs_nestor_params_path` is the absolute path to your `nestor_params.yaml` file.

**Note** One can also compare two completed NestOR runs with different parameter settings by running `python utils/compare_runs_v2_w_pyplot.py comparison_title run_set1 run_set2 ...` where comparison_title is the title for the runs to be compared, run_set1 and run_set2 are the NestOR runs to be compared.

## Outputs 

### Plots
Once the runs terminate, these plots will be generated: 
1. **Evidence**: The plot shows the mean values of evidence for all the candidate representations along with errorbars showing the standard error on the mean.
2. **MCMC per-step time**: The plot shows the time required to sample one MCMC step per run. This is computed as `(time taken for iteration 0)/((number of initial frames)*(number of MCMC steps per frame))`
3. **NestOR efficiency**: A scatter plot showing the mean total NestOR process time for all candidate representations.  

### Output YAML file

**Evidence related:**  
- log_estimated_evidence: `float`  
    _The estimated evidence value represented as natural logarithm of the estimated evidence_
- obtained_information: `float`  
    _Information obtained from the nested sampling run_
- analytical_uncertainty: `float`  
    _The analytical uncertainty associated with a evidence estimation for a run by nested sampling_

**Efficiency related**   
- mcmc_step_time: `float`  
    _Time taken per MCMC step. This is computed as `(time taken for iteration 0)/((number of initial frames)*(number of MCMC steps per frame))`_
- nestor_process_time: `float`  
    _Wall clock time taken by a nested sampling run to finish, represented in seconds_

**Termination related** 
- exit_code: `int` (0, 11, 12, 13)  
    _Exit code for a nested sampling run_
- termination_mode: `str`  
    _Cause for run termination_
- failed_iter: `int`  
    _Number of times Replica Exchange failed to obtain a sample from constrained prior in the current iteration of nested sampling_
- last_iter: `int`  
    _Iteration count (number of iterations) when nested sampling terminated_
- plateau_hits: `int`  
    _Number of consecutive times the nested sampling protocol detected a plateau in the estimated evidence_

**Exit codes:**  
- Exit code 0: Run terminated normally.  
- Exit code 11: Run terminated due to either a shuffle configuration error or NaN was encountered in the likelihoods. The run will be restarted automatically.  
- Exit code 12: Run terminated as NestOR ran out of maximum allowed iterations. The run will not be restarted.  
- Exit code 13: Run  terminated due to *Math domain error* in analytical uncertainty calculation. This happened probably because the run terminated too early resulting in a negative value for H. 

**_If a run terminates with `exit code = 12`, the run is considered incomplete (and is not rerun) and its results are not considered valid, i.e. these are not plotted and not used to infer optimal representation. Results from runs with exit codes 0 and 13 are used to infer the optimal representation_**

## Choice of NestOR parameters 

## **Information**
**Author(s):** Shreyas Arvindekar, Aditi Pathak, Kartik Majila, Shruthi Viswanath  
**Date**: April 7th, 2023  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Last known good IMP version:** `not tested`   
**Testable:** Yes  
**Parallelizeable:** Yes  
**Publications:**  Arvindekar, S, Pathak., A.S., Majila. K.M., Viswanath, S. Title to be decided. DOI: [](https://doi.org/).     
**_#TODO: Add title, DOI and link to the DOI_**
