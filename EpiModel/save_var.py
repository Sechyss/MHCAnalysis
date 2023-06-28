import numpy as np
import os


def save_vals(ResultPop, PathPop, UniqueAntigenCombos, UniqueHLAHaplos, Results, MHClookup_pre, initpop,
              foundingHLAhaplos, HLA_haplo_tbl, HLA_labels, SumLoc1rec, SumLoc2rec, Repeat, probMakeTcell):
    """
  This function saves the results of the model to a file.

  Args:
    ResultPop: The number of individuals in the population at each time step.
    PathPop: The number of infections at each time step.
    UniqueAntigenCombos: The number of deaths from infection at each time step.
    UniqueHLAHaplos: The number of individuals with antibodies to the pathogen at each time step.
    Results: The number of individuals with T cell responses to the pathogen at each time step.
    MHClookup_pre: The MHC lookup table.
    initpop: The population at the start of the simulation.
    foundingHLAhaplos: The founding HLA haplotypes.
    HLA_haplo_tbl: The HLA haplotype table.
    HLA_labels: The HLA labels.
    SumLoc1rec: The number of individuals with the first MHC locus.
    SumLoc2rec: The number of individuals with the second MHC locus.
    Repeat: The number of repeats of the simulation.
    probMakeTcell: The probability of an individual making a T cell response to the pathogen.

  Returns:
    None.
  """

    # Create the filename.

    filename = os.path.join(os.path.expanduser('~/OneDrive - University of Warwick/MatlabSimulations'),
                            str(probMakeTcell), 'repeat', str(Repeat), '.mat')

    # Save the results to the file.

    np.savez(filename, ResultPop=ResultPop, PathPop=PathPop, UniqueAntigenCombos=UniqueAntigenCombos,
             UniqueHLAHaplos=UniqueHLAHaplos, Results=Results, MHClookup_pre=MHClookup_pre,
             initpop=initpop, foundingHLAhaplos=foundingHLAhaplos, HLA_haplo_tbl=HLA_haplo_tbl,
             HLA_labels=HLA_labels, SumLoc1rec=SumLoc1rec, SumLoc2rec=SumLoc2rec)
