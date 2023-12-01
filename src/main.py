import argparse
import numpy as np

from chemicalio import importParameters, printInfo
from odetools import simulation
from errorsCheck import resetInfo

"""
parameters = allParameters[0]
#parameters = [chi, delta, ro, k, Da, As, div]
chi, delta, ro, k, Da, As, div = parameters

environment = allParameters [1]
# environment = [nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp]; 
nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = environment
"""

def printFinalInfo (parameters, environment, chemicalSpecies, reactions, matrixSimulation, timeSimulation): 

    print ("\n->END SIMULATION<-")
    printInfo(parameters, environment, chemicalSpecies, reactions)


def main(verbose, reset, file, importView):

    if reset: 
        resetInfo()

    parameters, environment, chemicalSpecies, reactions = importParameters (verbose, file)

    if importView: 
        printInfo(parameters, environment, chemicalSpecies, reactions)
        quit()

    if verbose:
        printInfo(parameters, environment, chemicalSpecies, reactions)

    (timeSimulation, matrixSimulation) = simulation (verbose, environment, parameters, chemicalSpecies, reactions)

    printFinalInfo (parameters, environment, chemicalSpecies, reactions, matrixSimulation, timeSimulation)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="ProtoGen Simulator v.1 - 284660")
    parser.add_argument("-v", "--verbose", action="store_true", help="additional prints")
    parser.add_argument("-r", "--reset", action="store_true", help="reset directory out/")
    parser.add_argument("-f", "--file", action="store_true", help="specify input file for parameters and reactions")
    parser.add_argument("-i", "--importV", action="store_true", help="view data imported")

    args = parser.parse_args()

    main(args.verbose, args.reset, args.file, args.importV)

