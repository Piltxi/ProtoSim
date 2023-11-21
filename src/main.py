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


def main(verbose, reset, file):

    if reset: 
        resetInfo()

    parameters, environment, chemicalSpecies, reactions = importParameters (verbose, file)

    if verbose:
        printInfo(parameters, environment, chemicalSpecies, reactions)

    (timeSimulation, matrixSimulation) = simulation (verbose, environment, parameters, chemicalSpecies, reactions)

    printFinalInfo (parameters, environment, chemicalSpecies, reactions, matrixSimulation, timeSimulation)

    excelExport(matrixSimulation, timeSimulation, chemicalSpecies, [parameters, environment, reactions])


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="ProtoGen Simulator v.1 - 284660")
    parser.add_argument("-v", "--verbose", action="store_true", help="additional prints")
    parser.add_argument("-r", "--reset", action="store_true", help="reset directory out/")
    parser.add_argument("-f", "--file", action="store_true", help="specify input file for parameters and reactions")

    args = parser.parse_args()

    main(args.verbose, args.reset, args.file)

""""
#! to do: 
1. copia dei parametri nei file di uscita
2. salvataggio dati su excel prima di chiusura
3. cambiare il sistema della memoria per memorizzare solo ultima riga di ultima gneerazione ( a scelta)
4. aggiunta flussi di controllo cstr con sintassi: 
>a; param (i)
b> ; param (k)

Si introducono due reazioni in più all'interno di ode, che simulano in autonomia il cstr. 
Sistemare generazione da espandere
"""