"""
#############################################
##  IMP Tutorial Script
##
#############################################
#
# Short modeling script combining EM and Crosslinking data
# to localize two domains of RNA Polymerase II
#
# Authors: Riccardo Pellarin, Charles Greenberg, Daniel Saltzberg
#
# References: Papers where this data is shown...
#
"""
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
# import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import ihm.cross_linkers

import sys,os


# ---------------------------
# Define Input Files
# ---------------------------
output_dir = 'run_'+sys.argv[1]
# flex_res = sys.argv[2]          # Resolution of flexible beads
init_frames = 10
nframes = 25 #int(sys.argv[2])

datadirectory = "../../../data/"
topology_file = datadirectory + sys.argv[2]
target_gmm_file = datadirectory+'emd_1883.map.mrc.gmm.50.txt'

# Initialize model
m = IMP.Model()

# Read in the topology file.
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                           pdb_dir=datadirectory,
                                           fasta_dir=datadirectory,
                                           gmm_dir=datadirectory)

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)

# Each state can be specified by a topology file.
bs.add_state(topology)

# Build the system representation and degrees of freedom
root_hier, dof = bs.execute_macro(max_rb_trans=4.0,
                                  max_rb_rot=0.3,
                                  max_bead_trans=4.0,
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)

# Fix all rigid bodies but not Rpb4 and Rpb7 (the stalk)
# First select and gather all particles to fix.
fixed_particles = []
for prot in ["Rpb1", "Rpb2", "Rpb3", "Rpb5", "Rpb6", "Rpb8", "Rpb9",
             "Rpb10", "Rpb11", "Rpb12"]:
    fixed_particles += IMP.atom.Selection(
        root_hier, molecule=prot).get_selected_particles()

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads, fixed_rbs = dof.disable_movers(fixed_particles,
                                            [IMP.core.RigidBodyMover,
                                             IMP.pmi.TransformMover])

# Randomize the initial configuration before sampling, of only the molecules
# we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50,
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)

outputobjects = []  # reporter objects (for stat files)
nester_restraints=[]

# -----------------------------------
# Define Scoring Function Components
# -----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated
#  into scoring. We also add them to the outputobjects list, so they
#  are reported in stat files.

# Connectivity keeps things connected along the backbone (ignores if inside
# same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname = mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
        mol, scale=2.0)
    cr.add_to_model()
    cr.set_label(molname)
    outputobjects.append(cr)

# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    included_objects=root_hier, resolution=10)
ev.add_to_model()
outputobjects.append(ev)

# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare
#  model to data. First the EM map is approximated with a Gaussian
#  Mixture Model (done separately). Second, the components of the model
#  are represented with Gaussians (forming the model GMM).
#  Other options:
#    scale_to_target_mass ensures the total mass of model and map are identical
#    slope: nudge model closer to map when far away
#    weight: experimental, needed becaues the EM restraint is quasi-Bayesian
em_components = IMP.pmi.tools.get_densities(root_hier)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                                 target_gmm_file,
                                                 scale_target_to_mass=True,
                                                 slope=0.000001,
                                                 weight=0)      #80.0)
gemt.add_to_model()
outputobjects.append(gemt)
nester_restraints.append(gemt)
# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging
#  components together)

# ---------
# Sampling XLs

xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_standard_keys()

xl1db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl1db.create_set_from_file("../../../data/sampling_xls.dat")

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xl1db,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="Trnka",
                                   linker=ihm.cross_linkers.dss,
                                   weight=1)

xl1.add_to_model()             # crosslink must be added to the model
outputobjects.append(xl1)
# nester_restraints.append(xl1)


# Crosslinks - dataset 2
#  We can easily add a second set of crosslinks.
#  These have a different format and label, but other settings are the same

# ---------
# Nesting XLs
xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_standard_keys()

xl2db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl2db.create_set_from_file("../../../data/nesting_xls.dat")

xl2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   database=xl2db,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="Chen",
                                   linker=ihm.cross_linkers.bs3,
                                   weight=0)
xl2.add_to_model()
outputobjects.append(xl2)
nester_restraints.append(xl2)

for i in outputobjects:
    print(i)
print("Nesting using these restraints:")
for j in nester_restraints:
    print(j)
# --------------------------
# Monte-Carlo Sampling
# --------------------------

# --------------------------
# Set MC Sampling Parameters
# --------------------------
# if '--test' in sys.argv:
#     num_frames = 100

# This object defines all components to be sampled as well as the
# sampling protocol
mc1 = IMP.pmi.macros.ReplicaExchange0(
    m, root_hier=root_hier, monte_carlo_sample_objects=dof.get_movers(),
    output_objects=outputobjects, monte_carlo_temperature=1.0,
    # simulated_annealing=True,
    # simulated_annealing_minimum_temperature=1.0,
    # simulated_annealing_maximum_temperature=2.5,
    # simulated_annealing_minimum_temperature_nframes=200,
    # simulated_annealing_maximum_temperature_nframes=20,
    replica_exchange_minimum_temperature=1.0,
    replica_exchange_maximum_temperature=2.5,
    number_of_best_scoring_models=0,
    monte_carlo_steps=10,
    number_of_frames=1,
    global_output_directory=output_dir,
    use_nester=True,
    nester_restraints=nester_restraints)

# Start Sampling
# mc1.execute_macro()
ns = IMP.pmi.macros.NestedSampling(num_init_frames=init_frames,
                                   num_frames_per_iter = nframes,
                                   nester_niter=5000,
                                   rex_macro=mc1,
                                   nester_restraints=nester_restraints,
                                   stopper_cap=20,
                                   early_stopper_cap=20)
ns.execute_nested_sampling()