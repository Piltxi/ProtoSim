import argparse

import chemicalio
from odetools import simulation

def main(verbose):

    parameters, environment, chemicalSpecies, reactions = chemicalio.importParameters (verbose)

    if verbose:
        chemicalio.printInfo(parameters, environment, chemicalSpecies, reactions)

    (mat_, time) = simulation (verbose, environment, parameters, chemicalSpecies, reactions)

    # chemicalio.excelExport(mat_, time_)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="ProtoGen Simulator v.1 - 284660")
    parser.add_argument("-v", "--verbose", action="store_true", help="additional prints")

    args = parser.parse_args()

    main(args.verbose)

