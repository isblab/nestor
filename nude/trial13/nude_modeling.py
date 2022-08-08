
##########################################################################################
######################### IMP Modeling Script for NuDe Complex ###########################
##########################################################################################


# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
#import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
#import IMP.saxs
import os
import sys
import pickle


# runType = 'test' # Specify test or prod
runID = sys.argv[1]   # Specify the number of runs
nframes = 25
# if runType == "test":
#     num_frames = 5000
# elif runType == "prod":
#     num_frames = 20000

max_shuffle_core = 5
max_shuffle_set2 = 50
rex_max_temp = 2.4

# Identify data files
sampling_adh_xl_data = "../../../input/xlms/sampling_adh_master.dat"
sampling_bs3dss_xl_data = "../../../input/xlms/sampling_bs3dss_master.dat"
sampling_dmtmm_xl_data =  "../../../input/xlms/sampling_dmtmm_master.dat"

evi_adh_xl_data = "../../../input/xlms/evicalc_adh_master.dat"
evi_bs3dss_xl_data = "../../../input/xlms/evicalc_bs3dss_master.dat"
evi_dmtmm_xl_data =  "../../../input/xlms/evicalc_dmtmm_master.dat"

gmm_data = "../../../input/gmm/emd_22904.gmm.40.txt"
dat_dir = "../../../input/"
# Restraint weights
intra_xl_weight = 1.0
inter_xl_weight = 10.0
xl_weight = 10
em_weight = 1000.0

# Topology File
topology_file = dat_dir + sys.argv[2]                                                      #"../../input/topology100.txt"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)

# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans= 1,
                                  max_rb_rot= 0.1,
                                  max_bead_trans= 3.2,
                                  max_srb_trans= 0.01,
                                  max_srb_rot=0.04)

################################################################################
########################## Fixing Particles ####################################
################################################################################
# # First select and gather all particles to fix.
# fixed_particles=[]
# for prot in ["MTA1"]:
#     for cp in [0,1]:
#         fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot,copy_index=cp,residue_indexes=range(165,334)).get_selected_particles()
# for prot in ["HDAC1"]:
#     for cp in [0,1]:
#         fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot,copy_index=cp,residue_indexes=range(8,377)).get_selected_particles()
#
# print(fixed_particles)
# # Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# # The flexible beads will still be flexible (fixed_beads is an empty list)!
# fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
#                                          [IMP.core.RigidBodyMover,
#                                           IMP.pmi.TransformMover])
# print('######################################')
# print(fixed_beads)


set1_core = []
for prot in ["MTA1"]:
    for cp in [0,1]:
        set1_core += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(165,334)).get_selected_particles()
for prot in ["HDAC1"]:
    for cp in [0,1]:
        set1_core += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(8,377)).get_selected_particles()
# print(set1_core)

set2 = []
for prot in ["MTA1"]:
    for cp in [0,1]:
        set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(1,165)).get_selected_particles()
for prot in ["MTA1"]:
    for cp in [0,1]:
        set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(334,432)).get_selected_particles()
for prot in ["HDAC1"]:
    for cp in [0,1]:
        set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(1,8)).get_selected_particles()
for prot in ["HDAC1"]:
    for cp in [0,1]:
        set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(377,483)).get_selected_particles()
for prot in ["MBD3"]:
    set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=0, residue_indexes=range(1,296)).get_selected_particles()
for prot in ["P66A"]:
    set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=0, residue_indexes=range(136,179)).get_selected_particles()
for prot in ["RBBP4"]:
    for cp in [0,1,2,3]:
        set2 += IMP.atom.Selection(root_hier, molecule=prot, copy_index=cp, residue_indexes=range(1,425)).get_selected_particles()
# print(set2)

fixed_set1_core_beads,fixed_set1_core = dof.disable_movers(set1_core,
                                                            [IMP.core.RigidBodyMover, IMP.pmi.TransformMover])

fixed_set2_beads, fixed_set2 = dof.disable_movers(set2,
                                 [IMP.core.RigidBodyMover, IMP.pmi.TransformMover])


# print(fixed_set1_core)
# print(fixed_set2)
# It's useful to have a list of the molecules.
molecules = t.get_components()


##### Uncomment the following lines to get test.rmf file to visualise the system representation
#
# # Uncomment this line for verbose output of the representation
# IMP.atom.show_with_representations(root_hier)
# # output to RMF
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)




#####################################################
##################### RESTRAINTS ####################
#####################################################

output_objects = []
nester_restraints = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print("Connectivity restraint applied")


# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT

evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)

print("Excluded volume restraint applied")

# -------------------------
# %%%%% CROSSLINKING RESTRAINT
# Sampling XL restraints

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xldb_adh_sampling = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_adh_sampling.create_set_from_file(file_name=sampling_adh_xl_data,
                              converter=xldbkc)
xlr_adh_sampling = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_adh_sampling, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,          # This adds a linear term to the scoring function
                label="adh",                        #   to bias crosslinks towards each other
                weight=xl_weight)       # Scaling factor for the restraint score.

xldb_bs3dss_sampling = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_bs3dss_sampling.create_set_from_file(file_name=sampling_bs3dss_xl_data,
                                 converter=xldbkc)
xlr_bs3dss_sampling = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_bs3dss_sampling, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="bs3dss",                        #   to bias crosslinks towards each other
                weight=xl_weight)       # Scaling factor for the restraint score.

xldb_dmtmm_sampling = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_dmtmm_sampling.create_set_from_file(file_name=sampling_dmtmm_xl_data,
                                converter=xldbkc)
xlr_dmtmm_sampling = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_dmtmm_sampling, # The crosslink database.
                length=16,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="dmtmm",                        #   to bias crosslinks towards each other
                weight=xl_weight)       # Scaling factor for the restraint score.


output_objects.append(xlr_adh_sampling)
output_objects.append(xlr_bs3dss_sampling)
output_objects.append(xlr_dmtmm_sampling)

# ----------------------------------------------------
# Nesting XL restraints

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xldb_adh_evicalc = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_adh_evicalc.create_set_from_file(file_name=evi_adh_xl_data,
                              converter=xldbkc)
xlr_adh_evicalc = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_adh_evicalc, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,          # This adds a linear term to the scoring function
                label="adh",                        #   to bias crosslinks towards each other
                weight=0)       # Scaling factor for the restraint score.

xldb_bs3dss_evicalc = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_bs3dss_evicalc.create_set_from_file(file_name=evi_bs3dss_xl_data,
                                 converter=xldbkc)
xlr_bs3dss_evicalc = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_bs3dss_evicalc, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="bs3dss",                        #   to bias crosslinks towards each other
                weight=0)       # Scaling factor for the restraint score.

xldb_dmtmm_evicalc = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_dmtmm_evicalc.create_set_from_file(file_name=evi_dmtmm_xl_data,
                                converter=xldbkc)
xlr_dmtmm_evicalc = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_dmtmm_evicalc, # The crosslink database.
                length=16,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="dmtmm",                        #   to bias crosslinks towards each other
                weight=0)       # Scaling factor for the restraint score.


output_objects.append(xlr_adh_evicalc)
output_objects.append(xlr_bs3dss_evicalc)
output_objects.append(xlr_dmtmm_evicalc)

nester_restraints.append(xlr_adh_evicalc)
nester_restraints.append(xlr_bs3dss_evicalc)
nester_restraints.append(xlr_dmtmm_evicalc)

print("Cross-linking restraint applied")

# # -------------------------
# %%%%% EM RESTRAINT

densities = IMP.atom.Selection(root_hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()

emr = IMP.pmi.restraints.em.GaussianEMRestraint(
             densities,                 # Evaluate the restraint using these model densities
             target_fn=gmm_data,        # The EM map, approximated as a gaussian mixture model (GMM)
             slope=0.00000001,          # a small linear restraint to pull objects towards the EM map center
             scale_target_to_mass=True, # Normalizes the total density of the model wrs: EM map. Only set to true
                                        #   if the EM map and "densities" contain the same objects.
             weight=0)          # the scaling factor for the EM score

output_objects.append(emr)
nester_restraints.append(emr)
print("EM Restraint Applied")

#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
# print("The type of run is: " + str(runType))
# print("Number of sampling frames: " + str(num_frames))
# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=max_shuffle_set2,
                                    excluded_rigid_bodies=fixed_set1_core)
                                    # hierarchies_included_in_collision=fixed_set1_core)

IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=max_shuffle_core,
                                    excluded_rigid_bodies=fixed_set2)
                                    # hierarchies_included_in_collision=fixed_set2)
# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
print('I crash after this')

dof.optimize_flexible_beads(15)#00)

IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)

# Now, add all of the other restraints to the scoring function to start sampling
evr.add_to_model()
emr.add_to_model()
xlr_adh_sampling.add_to_model()
xlr_bs3dss_sampling.add_to_model()
xlr_dmtmm_sampling.add_to_model()

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        crosslink_restraints=[xlr_adh_sampling,xlr_bs3dss_sampling,xlr_dmtmm_sampling],
        # This allows viewing the crosslinks in Chimera. Also, there is not inter-protein ADH crosslink available. Hence it is not mentioned in this list
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	    monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory="run_output_dir",      # The output directory for this sampling run.
        output_objects=output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames=nframes,#)            # Total number of frames to run / write to the RMF file.
        use_nester=True)#,
        # nester_restraints=nester_restraints)
        #test_mode=test_mode)                    # (Ignore this) Run in test mode (don't write anything)

# Ok, now we finally do the sampling!
# rex.execute_macro()

# Outputs are then analyzed in a separate analysis script.

ns = IMP.pmi.macros.NestedSampling(num_init_frames=10,
                                   num_frames_per_iter = nframes,
                                   nester_niter=5000,
                                   rex_macro=rex,
                                   nester_restraints=nester_restraints,
                                   stopper_cap=25,
                                   early_stopper_cap=20)

ns.execute_nested_sampling()
