from __future__ import print_function
import numpy as np
import IMP
import IMP.core
from IMP.pmi.tools import get_restraint_set

class nested_sampler():
    # check that isd is installed
    try:
        import IMP.isd
        isd_available = True
    except ImportError:
        isd_available = False


    def __init__(self, model):
        self.model = model
        self.li_fname = "likelihoods.dat"
        self.monte_carlo_sample_objects = None
        self.vars = {}
        self.nester_restraints = None
        self.nester_niter = None

    def parse_likelihoods(self):
        likelihoods = []
        with open(self.li_fname,'r') as Lif:
            for line in Lif:
                likelihoods.append(np.float(line.strip()))
        self.likelihoods=likelihoods
        return likelihoods

    def get_new_sample_with_constraint(self,worst_likelihood):
        # print(f"Worst likelihood: {worst_likelihood}")
        sampler_mc = IMP.pmi.shreyas_samplers.MonteCarlo(
                        self.model, self.monte_carlo_sample_objects,
                        self.vars["monte_carlo_temperature"])
        
        curr_li = 0
        early_stopper_limit = 5_000
        es = 0
        while es < early_stopper_limit:
            sampler_mc.optimize(1)
            Li = 1
            for restraint in self.nester_restraints:
                Li = Li*restraint.get_likelihood()
            curr_li=Li
            
            if curr_li > float(worst_likelihood):
               return curr_li
            else:
                es += 1
                if es%100==0:
                    print(f"es_count: {es}\tWorst likelihood: {worst_likelihood}\tattempted_li: {Li}")        
        print("Ealy stopper triggered")
        return None

