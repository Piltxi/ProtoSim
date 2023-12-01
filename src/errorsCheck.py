import os
import subprocess
# from reactions import ReactionType

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
            
            for reaction in data[1]:
                reactants = reaction["in"]
                products = reaction["out"]

                for specie in reactants + products:
                    
                    if reaction["type"] == ReactionType.FLOWIN: 
                        if products[0] not in data [0]: 
                            print (f"\nERROR 03 - invalid reaction: '{products[0]}' chemical species unknown ")
                            quit()
                        else: 
                            continue

                    if reaction["type"] == ReactionType.FLOWOUT: 
                        if reactants[0] not in data [0]: 
                            print (f"\nERROR 03 - invalid reaction: '{reactants[0]}' chemical species unknown")
                            quit ()
                        else: 
                            continue

                    if specie not in data[0]: 
                        try: 
                            float (specie)
                        except ValueError:
                            print (f"\nERROR 03 - invalid reaction: '{specie}' chemical species unknown")
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
            if data[0] < 1 or data[0] > data[1]:
                print (f"\nERROR 07 - loading generation indexes to expand\nIf you don't want to export any specific generation, type '-1' in the parameters file.")
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