Here is an example describing the use of NestOR to compare six different coarse-grained representations of the Nucleosome Deacetylase (NuDe) sub-complex of the Nucleosome Remodeling and Deacetylase (NuRD) complex. The `input` and `nude_modeling.py` are adapted from the [Integrative model of the NuRD subcomplexes](https://github.com/isblab/nurd) repository.

**NOTE:** _Prior to running the example, please download EM map from EMD22094, extract it, rename it as `emd_22904.mrc` and place it in `input/gmm/` directory_

## Description of the example files
The `examples/input` comprises of the fasta sequences (`examples/input/fasta`), PDB structures (`examples/input/pdb`), crosslinking mass spectrometry data (`examples/input/xlms`) and an EM map along with the corresponding GMM representation (`examples/input/gmm`). Each target crosslink file from `examples/input/xlms/original_xl_data` is split into `sampling_` and `evicalc_` files as described in the main module `README`.

In addition it also comprises the `topology{x}.txt` files that define the representation to be used for running the `nude_modeling.py` modeling script. The `{x}` in the topology file's name corresponds to the number of amino acid residues to be coarse-grained to a single flexible bead for regions that lack a previously characterized structure. In this example, we are comparing the 1, 5, 10, 20, 30 and 50 residues per bead coarse-grained representations of the regions with unknown structure of NuDe sub-complex. It is to be noted that the regions with known structure will be modeled as 1 and 10 residues per bead representation. Importantly, note that *any* representation can be compared as long there is a valid modeling file associated with it (representation can be mentioned via a topology file or manually). 

It also contains the `nestor_params_optrep.yaml` file which defines the NestOR parameters. Each parameter is described in the file itself. 

The `nude_modeling.py` script is also adapted from the [Integrative model of the NuRD subcomplexes](https://github.com/isblab/nurd) repository for use with NestOR.

The `nestor_output.yaml` contains an example output for the given setup. In addition to this file, NestOR also saves a model from each iteration. These models are not included here due to space constraints. It also generates plots visualizing the log(evidence) (mean and standard error on the mean) (`examples/trial_optrep_params_evidence_errorbarplot.png`), MCMC per step sampling time (`examples/trial_optrep_params_persteptime.png`), NestOR total process time (`examples/trial_optrep_params_proctime.png`) and per step MCMC sampling time and log(evidence) together for all candidate representations (`examples/sterr_evi_and_proctime.png`). 

## Running nested sampling on this example
For this example, the user only needs to run 
```
{path_to_IMP_installation} python nude_modeling.py
``` 
This command will run one nested sampling run. However, for use with the `wrapper_v{x}.py` the user will need to use modify this `nude_modeling.py` a bit to use command-line arguments. In the current form, the command line arguments have been hard coded in `nude_modeling.py`. The user will need to make these `sys.argv` in the correct order for the wrapper to work. Please see the comments in the `nude_modeling.py` for more details.