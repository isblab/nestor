Here is an example describing the use of NestOR to compare six different coarse-grained representations of the Nucleosome Deacetylase (NuDe) sub-complex of the Nucleosome Remodeling and Deacetylase (NuRD) complex. The `input` and `nude_modeling.py` are adapted from the [Integrative model of the NuRD subcomplexes](https://github.com/isblab/nurd) repository.

## Description of the example files
The `example/input` comprises of the fasta sequences (`example/input/fasta`), PDB structures (`example/input/pdb`), crosslinking mass spectrometry data (`example/input/xlms`) and an EM map along with the corresponding GMM representation (`example/input/gmm`). Each of the target crosslink file from `example/input/xlms/original_xl_data` is split into `sampling_` and `evicalc_` files as described below.

In addition it also comprises the `topology{x}.txt` files that define the representation to be used for running the `nude_modeling.py` modeling script. The `{x}` in the topology file's name corresponds to the number of amino acid residues to be coarse-grained to a single flexible bead for regions that lack a previously characterized structure. In this example, we are comparing the 1, 5, 10, 20, 30 and 50 residues per bead coarse-grained representations of the regions with unknown structure of NuDe sub-complex. It is to be noted that the regions with known structure will be modeled as 1 and 10 residues per bead representation.

The `nude_modeling.py` script is also adapted from the [Integrative model of the NuRD subcomplexes](https://github.com/isblab/nurd) repository to implement NestOR.

The `nestor_params_optrep.yaml` defines the NestOR parameters. Each of the parameter is described in the file itself. 

The `nestor_output.yaml` contains an example output for the given setup. In addition to this file, NestOR also saves a model from each iteration. These models are not included here due to space constraints. It also generate teh plots visualizing the log(evidence) (mean and standard error on the mean) (`example/trial_optrep_params_evidence_errorbarplot.png`), MCMC per step sampling time (`example/trial_optrep_params_persteptime.png`) and NestOR total process time (`example/trial_optrep_params_proctime.png`) for all candidate representations. An additional script `figure_scripts/plot_evidence_proctime_together.py` can then be used to plot the per step MCMC sampling time and log(evidence together). The result from this will look similar to this: ![sterr_evi_and_proctime.png](https://github.com/isblab/nestor/blob/main/example/sterr_evi_and_proctime.png)

## Tutorial
1. Once the inputs, modeling script, and NestOR parameter file are made, run the NestOR wrapper as follows:
```python wrapper_v5.py {nestor_param_path} {mode}```
where, `nestor_param_path` refers to the absolute path to the `nestor_params.yaml` file and `mode` refers to the mode of representation (`manual`/`topology`). The default choice is `topology` (does not need to be mentioned), If the representation is defined in the modeling script, use `manual` argument.

2. Once all the runs are complete, run the following command:
```python figure_scripts/plot_evidence_proctime_together.py {path} {name}```
where, `path` refers to the `parent_dir` in NestOR params file and `name` refers to the name of the assembly for which the representations is being optimized.