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

### Outputs 

### Plots
Once the runs terminate, two plots will be generated: 
1. **Evidence**: The plot showing the mean values of evidence for all the candidate representations along with errorbars showing the standard error on the mean.
2. **MCMC per-step time**: #TODO
3. **NestOR efficiency**: A scatter plot showing the mean total NestOR process time for all candidate representations.  

### Output YAML file
#TODO describe the output values here. 

**Exit codes:**
**Exit code 0:** Run terminated normally.  
**Exit code 11:** Run terminated due to either a shuffle configuration error or NaN was encountered in the likelihoods. The run will be restarted automatically.  
**Exit code 12:** Run terminated as NestOR ran out of maximum allowed iterations. The run will not be restarted.  
**Exit code 13:** Run  terminated due to *Math domain error* in analytical uncertainty calculation. This happened probably because the run terminated too early resulting in a negative value for H. 

#TODO exit codes 12, 13: do we use in our statistics? Mention in README

## Choice of NestOR parameters 
#TODO Why make it redundant with the paper? 

## **Information**
**Author(s):** Shreyas Arvindekar, Aditi Pathak, Kartik Majila, Shruthi Viswanath  
**Date**: April 7th, 2023  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Last known good IMP version:** [![build info](https://integrativemodeling.org/systems/41/badge.svg?branch=main)](https://integrativemodeling.org/systems/)    #TODO: Fix this. Make it `not tested`   
**Testable:** Yes  
**Parallelizeable:** Yes  
**Publications:**  Arvindekar, S, Pathak., A.S., Majila. K.M., Viswanath, S. Title to be decided. DOI: [](https://doi.org/).     
**_#TODO: Add title, DOI and link to the DOI_**
