[![PubMed](https://salilab.org/imp-systems/static/images/pubmed.png)](https://pubmed.ncbi.nlm.nih.gov/xyz/)     #TODO: Replace xyz with PubMed ID


# **NestOR: Optimizing the representation for Integrative Structure Modeling using Nested Sampling**


### **Installation:**
#### **Dependencies:**  
Standard IMP installation (compiled from the source code), numpy, mergedeep, mpi4py, matplotlib, pyyaml, plotly
#### **NestOR installation:**  
1. Compile IMP from the souce code to your choice of directory
2. Replace the `macros.py` in `imp/modules/pmi/pyext/src/` with the macros.py in the current repository. Make sure the file is named `macros.py` in the destination directory. 
3. Similarly, replace the restraints directory in `imp/modules/pmi/pyext/src/` with the restraints directory in the present repository.


### **Running NestOR:**
1. Make the modeling script in the form as shown in the `example/modeling.py`. One will also need to make separate topology files for different candidate representations.
2. Set appropriate parameters in the `nestor_params.yaml` file. An example param file can be found in `examples/`.
3. Finally, run nestor with the following command `python wrapper_vX.py $abs_nestor_params_path` where wrapper_vX.py is the wrapper script in this repository and abs_nestor_params_path is the absolute path to your `nestor_params.yaml` file.

Once the runs terminate, two plots will be generated: 
1. The plot showing the mean values of evidence for all the candidate representations along with errorbars showing the standard error on the mean.
2. A scatter plot showing the mean total nestor process time for all candidate representations.  

One can also compare NestOR runs with different parameter settings by running `python utils/compare_runs_v2_w_pyplot.py comparison_title run_set1 run_set2 ...` where comparison_title is the title for the runs to be compared, run_set1 and run_set2 are tbe NestOR runs to be compared.


### **Exit codes:**
**Exit code 0:** Process terminated normally.  
**Exit code 11:** Process terminated due to either a shuffle configuration error or NaN was encountered in the likelihoods. The process will be restarted automatically.  
**Exit code 12:** Process terminated as NestOR ran out of maximum allowed iterations. The process will not be restarted.  
**Exit code 13:** Process terminated due to *Math domain error* in analytical uncertainty calculation. This happened probably because the run terminated too early resulting in a negative value for H.   


### **Information**
**Author(s):** Shreyas Arvindekar, Aditi Pathak, Kartik Majila, Shruthi Viswanath  
**Date**: April 7th, 2023  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Last known good IMP version:** [![build info](https://integrativemodeling.org/systems/41/badge.svg?branch=main)](https://integrativemodeling.org/systems/)    #TODO: Fix this. Make it `not tested`   
**Testable:** Yes  
**Parallelizeable:** Yes  
**Publications:**  Arvindekar, S, Jackman, MJ, Low, JKK, Landsberg, MJ, Mackay, JP, Viswanath, S. Title to be decided. DOI: [](https://doi.org/).     
**_TODO: Add title, DOI and link to the DOI_**