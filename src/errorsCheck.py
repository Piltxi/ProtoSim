import os
import subprocess

def _errorState (code, message): 
    print (f"\nERROR {code} - {message}")
    quit()

def checkProtoSim (arg, data):
    
    match arg: 

        case 0: 
            if (data != 'C'):
                print ("\nERROR 00 - loading parameters from file")
                quit()
        
        case 1: 
            print ("\nERROR 01 - information symbol: view chemistry file", data)
            quit ()
        
        case 2: 
            print (f"\nERROR 02 - \" {data} \" type of reaction unknown")
            quit()

        case 3: 
            
            from reactions import ReactionType
            from chemicalio import getOrdinal

            for i, reaction in enumerate (data[1]):
                
                reactants = reaction["in"]
                products = reaction["out"]

                for specie in reactants + products:
                    
                    if reaction["type"] == ReactionType.FLOWIN: 
                        if products[0] not in data [0]: 
                            print (f"\nERROR 03 - invalid {getOrdinal(i+1)} reaction: '{products[0]}' chemical species unknown ")
                            quit()
                        else: 
                            continue

                    if reaction["type"] == ReactionType.FLOWOUT: 
                        if reactants[0] not in data [0]: 
                            print (f"\nERROR 03 - invalid {getOrdinal(i+1)} reaction: '{reactants[0]}' chemical species unknown")
                            quit ()
                        else: 
                            continue

                    if reaction["type"] == ReactionType.DIFFUSION: 
                        
                        if products[0] not in data [0]: 
                            print (f"\nERROR 03 - invalid {getOrdinal(i+1)} reaction: '{products[0]}' chemical species unknown")
                            quit ()
                        else: 
                            continue

                    if specie not in data[0]: 
                            print (f"\nERROR 03 - invalid {getOrdinal(i+1)} reaction: '{specie}' chemical species unknown")
                            quit()

        case 4: 
            print ("\nERROR 04 - negative quantities detected\n")
            print ("Value: ", data[0], "\tIndex: ", data[1], "\n")
            quit()

        case 5: 
            print ("\nERROR 05 - ode_function, unknow variation rules for ", data["type"])
            quit()

        case 6: 
            print (f"\nERROR 06 - indexing reactions '{data}'")
            quit()

        case 7: 

            nIterates, gen_exp, genExp_timing = data

            """
                check list: 
                1] nIterates > 1
                2] expand gen not enabled
                3] expand index(s) consistency with nIterates
                4] expand index(s) consistency with timing parameter(s)
                5] timing parameter(s) > 0
                6] badref for expand index(s) enabled
            """

            # 1-> numbers of nIterates imported: <1 
            if nIterates < 1: 
                _errorState ([7, 0], "invalid number of iterations [nIterates] detected")

            # 2-> no expand request
            if len(gen_exp) == 1 and gen_exp[0] == -1: 
                if len(genExp_timing) == 1 and genExp_timing[0] == -1:
                    return
                else: 
                    _errorState ([7, 2], f"expansion indices and corresponding timing\ninconsistent values ​​detected between timing({genExp_timing}), {len(genExp_timing)}value(s) and expansions({gen_exp}), {len(gen_exp)}value(s)")
            
            # 3-> expansion indices corrected with nIterates
            for element in gen_exp: 
                if element <= 0 or element > nIterates:
                    _errorState ([7, 1], "loading generation indexes to expand\nIf you don't want to export any specific generation, type '-1' in the parameters file.")

            # 4-> expansion indices == timing indices
            if len(gen_exp) != len(genExp_timing): 
                _errorState ([7, 2], f"expansion indices and corresponding timing\ninconsistent values ​​detected between timing({genExp_timing}), {len(genExp_timing)}value(s) and expansions({gen_exp}), {len(gen_exp)}value(s)")

            # 5-> timing indices > 0
            i=1
            for element in genExp_timing: 
                if element <= 0 and not element == -1 :
                    _errorState ([7, 3], f"unknown time of expand n. {i}\nSpecify timing parameter for every generation to expand, or type '-1'")
                i+=1

            # 6-> badref for expand index(s) enabled
            if gen_exp[0] == -1 and genExp_timing[0] != -1:
                _errorState ([7, 4], f"unknown time of export generations of expansion\nIf you don't want to export any specific generation, type '-1' in the parameters file.")

            # for element in gen_exp: 
            #     element-=1

        case 8: 

            if data [0] <= 0:
                print (f"\nERROR 08 - zero reactions found\n")
                quit()

            if data [1] < 0 or not (isinstance(data [1], int)):

                print (f"\nERROR 08 - unknown flux number: please type '0' in parameters file to disable tracking\n")
                quit()

            if data [1] == 0:
                return

            if data [1] > 0: 
                if data [1] > data [0]: 
                    print (f"\nERROR 08 - too many fluxes recognized\nnumber of imported reactions: {data [0]} - number of flux imported: {data[1]}.")
                    quit()

        case 9: 
            print (f"\nERROR 09 - tollerance Test: empty protoX [{data}]\n")
            quit()

        case 10:
            print (f"\nERROR 10 - species '{data}' not found in loadedSpecies\n")
            quit()

        case _: 
            print ("\nUNKNOW ERROR XY")
            quit ()

def resetInfo (): 
    
    if os.path.exists("../out"):
        try:
            subprocess.run(["rm", "-fr", "../out"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error in removing the existing directory: {e}")
    quit()