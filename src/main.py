import argparse
import numpy as np

from chemicalio import importParameters, excelExport, printInfo
from odetools import simulation
from errorsCheck import resetInfo

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
    printInfo(parameters, environment, chemicalSpecies, reactions)


def main(verbose, reset):

    if reset: 
        resetInfo()

    parameters, environment, chemicalSpecies, reactions = importParameters (verbose)

    if verbose:
        printInfo(parameters, environment, chemicalSpecies, reactions)

    (timeSimulation, matrixSimulation) = simulation (verbose, environment, parameters, chemicalSpecies, reactions)

   # printFinalInfo (parameters, environment, chemicalSpecies, reactions, matrixSimulation, timeSimulation)

    excelExport(matrixSimulation, timeSimulation, chemicalSpecies, [parameters, environment])


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="ProtoGen Simulator v.1 - 284660")
    parser.add_argument("-v", "--verbose", action="store_true", help="additional prints")
    parser.add_argument("-r", "--reset", action="store_true", help="reset out/ directory")

    args = parser.parse_args()

    main(args.verbose, args.reset)

