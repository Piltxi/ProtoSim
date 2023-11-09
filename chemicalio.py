import numpy as np
import subprocess
import os

from reactions import ReactionType, identifyType
import errorsCheck

def importParameters (verbose): 
    
    fi=open("parameters.txt",'r')   

    delta = eval(fi.readline().split()[0])
    ro = eval(fi.readline().split()[0])
    Da = eval (fi.readline().split()[0])
    k = eval(fi.readline().split()[0])
    As = eval(fi.readline().split()[0])
    div = eval(fi.readline().split()[0]) 
    nIterates = eval(fi.readline().split()[0])
    t_end = eval(fi.readline().split()[0])
    max_step = eval(fi.readline().split()[0])  
    toll_min = eval(fi.readline().split()[0])
    toll_max = eval(fi.readline().split()[0])
    nFlux = eval(fi.readline().split()[0]) 
    nGen = eval(fi.readline().split()[0])

    fi.close()

    calving = 0.353553
    chi = 1/(6*pow(np.pi*pow(delta,3)*pow(ro,3),0.5))
    # print ("chi main: ", chi)
    t_start = 0.

    # List of parameters to resolve ODE
    parameters = [chi, delta, ro, k, Da, As, div]

    # List of environment sets
    environment = [nIterates, t_end, max_step, toll_min, toll_max, nFlux, nGen, calving]; 
    
    chemicalSpecies = {}
    reactions = []

    fi=open("chimica.txt",'r')
    specify = fi.readlines()
    fi.close()

    for line in specify: 
        
        line = line.strip()
        
        if line: 
            parts = line.split ("\t")
            
            if len(parts) == 3: # if  parts.strip (): if ';' in parts: 
                species, quantity, coefficient = parts

                if  verbose:
                    print (f"Loaded chemical species: {species}\t{quantity} Kg\t{coefficient}")

                chemicalSpecies[species] = (float(quantity), float(coefficient))
            
            else:  
                
                if verbose:
                    print ("start import", end = "|\t") 
                
                reactionsParts = line.strip().split(";")
                
                if len (reactionsParts) == 2: 
                    
                    reactionType = ReactionType.ND

                    reaction_str = reactionsParts[0].strip()
                    coefficient = float(reactionsParts[1].strip())
                    reagents, products = reaction_str.split('>')
                    reagents = [reagent.strip() for reagent in reagents.split('+')]
                    products = [product.strip() for product in products.split('+')]
                    
                    reactionType = identifyType(line, verbose)

                    reaction_data = {
                        "in": reagents,
                        "out": products,
                        "k": coefficient,
                        "type": reactionType
                    }
                    
                    if verbose: 
                        print ("reactants: ", reagents, "|products: ", products, "|type: ", reactionType.name, "|Coefficient: ", coefficient, end=" \t|")
                        print ("end import")
                    
                    reactions.append(reaction_data)
                
                else: 
                    checkProtoSim (1, "error")
    
    loadedSpecies = list(chemicalSpecies.keys())
    
    # Check of chemical species in the reactions
    errorsCheck.checkProtoSim(3, [loadedSpecies, reactions])
    
    
    return [parameters, environment, chemicalSpecies, reactions]

def printInfo (parameters, environment, chemicalSpecies, reactions): 

    sAlpha = "\u03B1"
    sEta = "\u03b7"
    sDelta = "\u03b4"
    sMu = "\u03bc"
    sRo = "\u03c1"
    sChi = "\u03c7"
    sPi = "\u03c0"
    sNu = "\u03bd"
    sBChi = "\u03a7"
    sKappa = "\u039a"

    chi, delta, ro, k, Da, As, div = parameters
    nIterates, t_end, max_step, toll_min, toll_max, nFlux, nGen, calving = environment

    print ("\nSTART PRINTING INFO...\n")

    print("Recognized Parameters:")
    print(sDelta, ":\t", delta)
    print(sRo, ":\t", ro)
    print (sChi, ":\t", chi)
    print ("As:\t", As)
    print (sKappa,":\t", k)
    print("Da:\t", Da)

    print("\nExecution Parameters:")
    print("flux:\t",nFlux)
    print("generations:\t",nGen)
    print("iterations:\t", nIterates)
    print("calving:\t", calving)
    print("end time:\t",t_end)
    print("max  step:\t",max_step)       
    print("min toll. :",toll_min, "\tmax toll. :",toll_max)

    print("\nChemical Species imported:")
    i = 0
    for species, (quantity, coefficient) in chemicalSpecies.items():
        i+=1
        print (f"{i}] {species} \t {quantity} kg\tCoefficient: ", coefficient)

    print("\nReactions imported:")
    i = 0
    for i, reaction in enumerate (reactions): 
        reagents_str = " + ".join(reaction["in"])
        products_str = " + ".join(reaction["out"])
        print (f"{i+1}] {reagents_str} -> {products_str}\nKinetic Coefficient: ", reaction["k"], "\tType: ", reaction["type"].value, "\n")

    protoGen = np.array([chemicalSpecies[quantity][0] for quantity in chemicalSpecies])
    loadedSpecies = list(chemicalSpecies.keys())
    
    # mapReactions = map_species_to_indices(reactions, protoGen, loadedSpecies)
    # print("\nReactions [index print]:")
    # printMapReactions(mapReactions)

    print ("\nEND PRINTING INFO...\n")

def map_species_to_indices(reactions, protoGen, loadedSpecies):
    index_based_reactions = []

    for reaction in reactions:
        indexed_reagents = []
        indexed_products = []

        for reagent in reaction["in"]:
            if reaction["type"] == ReactionType.DIFFUSION:
                indexed_reagents.append(float(reaction["in"][0]))
            elif reagent in loadedSpecies:
                indexed_reagents.append(loadedSpecies.index(reagent))

        for product in reaction["out"]:
            if product in loadedSpecies:
                indexed_products.append(loadedSpecies.index(product))
            else:
                checkProtoSim (6, product)

        indexed_reaction = {
            "in": indexed_reagents,
            "out": indexed_products,
            "k": reaction["k"],
            "type": reaction["type"]
        }

        index_based_reactions.append(indexed_reaction)
    
    return index_based_reactions

def printMapReactions (mapReactions): 
    for reaction in mapReactions:
        print("Reagents (indices):", reaction["in"])
        print("Products (indices):", reaction["out"])
        print("Rate constant (k):", reaction["k"])
        print("Reaction type:", reaction["type"])
        print()

def excelExport (mat, time): 

    # Sostituisci "nome_della_tua_nuova_directory" con il nome desiderato per la nuova directory
    directory_name = "out"

    # Verifica se la directory esiste già
    if os.path.exists(directory_name):
        # Rimuovi la directory esistente
        try:
            subprocess.run(["rm", "-r", directory_name], check=True, shell=True)
            print(f"La directory '{directory_name}' esistente è stata rimossa.")
        except subprocess.CalledProcessError as e:
            print(f"Errore durante la rimozione della directory esistente: {e}")

    # Ora puoi creare la nuova directory
    try:
        subprocess.run(["mkdir", directory_name], check=True, shell=True)
        print(f"La directory '{directory_name}' è stata creata con successo.")
    except subprocess.CalledProcessError as e:
        print(f"Errore durante la creazione della directory: {e}")

    name = f"out/out.xlsx"
    workbook = xlsxwriter.Workbook (name)
    wq = workbook.add_worksheet("Quantity")
    wc = workbook.add_worksheet("Concentration")

