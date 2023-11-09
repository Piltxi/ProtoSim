import os
import subprocess

def checkProtoSim (arg, data):
    
    match arg: 

        case 0: 
            if (data != 'C'):
                print ("\nERROR 00 - loading parameters from file")
                quit()
        
        case 1: 
            print ("\nERROR 01 - chemical symbol in reactions")
            quit ()
        
        case 2: 
            print (f"\nERROR 02 - \" {data} \" type of reaction unknown")
            quit()

        case 3: 

            for reaction in data[1]:
                reactants = reaction["in"]
                products = reaction["out"]
            
                for specie in reactants + products:
                    if specie not in data[0]: 
                        try: 
                            float (specie)
                        except ValueError:
                            print ("\nERROR 03 - invalid reaction: chemical species unknown")
                            quit()

        case 4: 
            print ("\nERROR 04 - negative quantities detected")
            print ("Data: ", data)
            quit()

        case 5: 
            print ("\nERROR 05 - ode_function error")
            print ("Data: ", data)

        case 6: 
            print ("\nERROR 06 - conversion dictionary(product unknow)")
            print ("Data: ", data)

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