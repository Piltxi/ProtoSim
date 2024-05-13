import os
import subprocess

def _errorState (code, message): 
    print (f"\nERROR {code} - {message}")
    quit()

def checkProtoSim (arg, data):
    
    match arg: 

        # potentially deprecated
        case 0: 
            if (data != 'C'):
                _errorState (0, "loading parameters from file")
        
        case 1: 
            _errorState (1, f"information symbol |{data}| loaded from chemistry file")
        
        case 2: 
            _errorState (2, f"|{data}| type of reaction unknown")

        case 3: 
            
            from reactions import ReactionType
            from chemicalio import getOrdinal

            for i, reaction in enumerate (data[1]):
                
                reactants = reaction["in"]
                products = reaction["out"]

                for specie in reactants + products:
                    
                    if reaction["type"] == ReactionType.FLOWIN: 
                        if products[0] not in data [0]: 
                            _errorState(3, f"invalid {getOrdinal(i+1)} reaction: '{products[0]}' chemical species unknown")
                        else: 
                            continue

                    if reaction["type"] == ReactionType.FLOWOUT: 
                        if reactants[0] not in data [0]: 
                            _errorState(3, f"invalid {getOrdinal(i+1)} reaction: '{reactants[0]}' chemical species unknown")
                        else: 
                            continue

                    if reaction["type"] == ReactionType.DIFFUSION: 
                        
                        if products[0] not in data [0]: 
                            _errorState(3, f"invalid {getOrdinal(i+1)} reaction: '{products[0]}' chemical species unknown")
                        else: 
                            continue

                    if specie not in data[0]: 
                            _errorState(3, f"invalid {getOrdinal(i+1)} reaction: '{specie}' chemical species unknown")

        case 4: 
            _errorState(4, "negative quantities detected\nvalue: {data[0]}\tindex: {data[1]}\n")

        case 5: 
            _errorState(5, f"unknow ODE rules for {data["type"]}")

        case 6: 
            _errorState(6, f"indexing reactions |{data}|")

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
                    _errorState ([7, 1], f"expansion indices and corresponding timing\ninconsistent values ​​detected between timing({genExp_timing}), {len(genExp_timing)}value(s) and expansions({gen_exp}), {len(gen_exp)}value(s)")
            
            # 3-> expansion indices corrected with nIterates
            for element in gen_exp: 
                if element <= 0 or element > nIterates:
                    _errorState ([7, 2], "loading generation indexes to expand\nIf you don't want to export any specific generation, type '-1' in the parameters file.")

            # 4-> expansion indices == timing indices
            if len(gen_exp) != len(genExp_timing): 
                _errorState ([7, 3], f"expansion indices and corresponding timing\ninconsistent values ​​detected between timing({genExp_timing}), {len(genExp_timing)}value(s) and expansions({gen_exp}), {len(gen_exp)}value(s)")

            # 5-> timing indices > 0
            i=1
            for element in genExp_timing: 
                if element <= 0 and not element == -1 :
                    _errorState ([7, 4], f"unknown time of expand n. {i}\nSpecify timing parameter for every generation to expand, or type '-1'")
                i+=1

            # 6-> badref for expand index(s) enabled
            if gen_exp[0] == -1 and genExp_timing[0] != -1:
                _errorState ([7, 5], f"unknown time of export generations of expansion\nIf you don't want to export any specific generation, type '-1' in the parameters file.")

            # for element in gen_exp: 
            #     element-=1

        case 8: 

            if data [0] <= 0:
                _errorState (8, "zero reactions found")

            if data [1] < 0 or not (isinstance(data [1], int)):
                _errorState (8, "unknown flux number: please type '0' in parameters file to disable tracking")

            if data [1] == 0:
                return

            if data [1] > 0: 
                if data [1] > data [0]: 
                    _errorState (8, f"too many fluxes recognized:\n imported reactions: |{data [0]}|\timported flux: |{data[1]}|")

        case 9: 
            _errorState (9, f"tollerance test: empty protoX {data}")

        case 10:
            _errorState (10, f"species |{data}| not found in loadedSpecie")

        case 11:
           # recived data = [gen_exp, genExp_time, target])
           _errorState (11, f"generation expansion: timing match not found\nGen. to expand: {data[0]}\nTiming: {data[1]}\nTarget: {data[2]}\n")

        case _: 
            _errorState ("XY", "UNKNOW ERROR")

def resetInfo (): 
    
    if os.path.exists("../out"):
        try:
            subprocess.run(["rm", "-fr", "../out"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error in removing the existing directory: {e}")
    quit()