import argparse
import numpy as np

import chemicalio
from odetools import simulation

"""
parameters = allParameters[0]
#parameters = [chi, delta, ro, k, Da, As, div]
chi, delta, ro, k, Da, As, div = parameters

environment = allParameters [1]
# environment = [nIterates, t_end, max_step, toll_min, toll_max, nFlux, nGen]; 
nIterates, t_end, max_step, toll_min, toll_max, nFlux, nGen, calving = environment
"""

def printFinalInfo (parameters, environment, chemicalSpecies, reactions, matrixSimulation, timeSimulation): 

    print ("End of simulation...\n")
    chemicalio.printInfo(parameters, environment, chemicalSpecies, reactions)





def main(verbose):

    parameters, environment, chemicalSpecies, reactions = chemicalio.importParameters (verbose)

    if verbose:
        chemicalio.printInfo(parameters, environment, chemicalSpecies, reactions)

    (matrixSimulation, timeSimulation) = simulation (verbose, environment, parameters, chemicalSpecies, reactions)
    printFinalInfo (parameters, environment, chemicalSpecies, reactions, matrixSimulation, timeSimulation)


    # chemicalio.excelExport(mat_, time_)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="ProtoGen Simulator v.1 - 284660")
    parser.add_argument("-v", "--verbose", action="store_true", help="additional prints")

    args = parser.parse_args()

    main(args.verbose)

