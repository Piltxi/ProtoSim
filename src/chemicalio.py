import numpy as np
import subprocess
import os
import xlsxwriter
from datetime import datetime

from errorsCheck import checkProtoSim
from reactions import identifyType, ReactionType

"""
parameters = allParameters[0]
#parameters = [chi, delta, ro, Da, div]
chi, delta, ro, Da, div = parameters

environment = allParameters [1]
# environment = [nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp]; 
nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = environment
"""

def getConcentration(x, C, ro, delta):
    return x / getVolume(C, ro, delta)

def getVolume(C, ro, delta):
    radius = delta * (pow(C / (ro * np.pi * delta * delta * delta) - 1/3, 0.5) - 1) / 2
    return 4 * np.pi * radius * radius * radius / 3

def importParameters (verbose, file): 

    parametersFile = "../input/"
    reactionsFile = "../input/"
    
    if file: 
        parametersUser = input("Type the name of text file with parameters> ")
        reactionsUser = input("Type the name of text file with reaction> ")
        parametersFile = parametersFile + parametersUser
        reactionsFile = reactionsFile + reactionsUser
    else: 
        parametersFile = parametersFile + "parameters.txt"
        reactionsFile = reactionsFile + "chimica.txt"

    fi=open(parametersFile,'r')   

    delta = eval(fi.readline().split()[0])
    ro = eval(fi.readline().split()[0])
    Da = eval (fi.readline().split()[0])
    div = eval(fi.readline().split()[0]) 
    nIterates = eval(fi.readline().split()[0])
    t_end = eval(fi.readline().split()[0])
    max_step = eval(fi.readline().split()[0])  
    toll_min = eval(fi.readline().split()[0])
    toll_max = eval(fi.readline().split()[0])
    nFlux = eval(fi.readline().split()[0]) 
    gen_exp = eval(fi.readline().split()[0])

    fi.close()

    checkProtoSim(7, [gen_exp, nIterates])
    gen_exp -= 1

    calving = 0.353553
    chi = 1/(6*pow(np.pi*pow(delta,3)*pow(ro,3),0.5))

    # List of parameters to resolve ODE
    parameters = [chi, delta, ro, Da, div]

    # List of environment sets
    environment = [nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving]; 
    
    chemicalSpecies = {}
    reactions = []

    fi=open(reactionsFile,'r')
    specify = fi.readlines()
    fi.close()

    for line in specify: 
        
        line = line.strip()
        
        if line: 
            #parts = line.split ("\t")
            
            parts = line.split()
            
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
                        print ("reactants: ", reagents, "|products: ", products, "|type: ", reactionType.value, "|Coefficient: ", coefficient, end=" \t|")
                        print ("end import")
                    
                    reactions.append(reaction_data)
                
                else: 
                    checkProtoSim (1, line)
    
    loadedSpecies = list(chemicalSpecies.keys())
     
    # Check of chemical species in the reactions
    checkProtoSim(3, [loadedSpecies, reactions])

    # Check acceptability of number of fluxes and number of reactions
    checkProtoSim (8, [len(reactions), nFlux])
    
    return [parameters, environment, chemicalSpecies, reactions]

def printInfo (parameters, environment, chemicalSpecies, reactions): 

    sDelta = "\u03b4"
    sRo = "\u03c1"
    sChi = "\u03c7"

    chi, delta, ro, Da, div = parameters
    nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = environment

    print("\nRecognized Parameters:")
    print(sDelta, ":\t", delta)
    print(sRo, ":\t", ro)
    print (sChi, ":\t", chi)
    print("Da:\t", Da)
    print("Duplication threshold:\t", div)

    print("\nExecution Parameters:")
    print("flux:\t",nFlux)
    print("iterations:\t", nIterates)
    print("calving:\t", calving)
    print("end time:\t",t_end)
    print("max  step:\t",max_step)       
    print("min toll. :",toll_min, "\tmax toll. :",toll_max)
    print("generation to expand:\t", gen_exp:=gen_exp+1)

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

        if reaction["type"] is not ReactionType.FLOWIN and reaction["type"] is not ReactionType.FLOWOUT:
            print (f"{i+1}] {reagents_str} -> {products_str}\nKinetic Coefficient: ", reaction["k"], "\tType: ", reaction["type"].value, "\n")

        else: 
            arrow =" \u2192 "
            if reaction["type"] is ReactionType.FLOWIN: 
                print (f"{i+1}] {products_str} {arrow} [CSTR] \nSubstance Fraction: ", reaction["k"], "\tType: ", reaction["type"].value, "\n")
            else: 
                print (f"{i+1}] [CSTR] {arrow} {reagents_str}\nSubstance Fraction: ", reaction["k"], "\tType: ", reaction["type"].value, "\n")

def map_species_to_indices(reactions, loadedSpecies):
    index_based_reactions = []

    for reaction in reactions:
        indexed_reagents = []
        indexed_products = []

        if reaction["type"] == ReactionType.FLOWIN or reaction["type"] == ReactionType.FLOWOUT: 

            if reaction["type"] == ReactionType.FLOWIN: 
                indexed_reagents.append(None)
                indexed_products.append(loadedSpecies.index(reaction["out"][0]))

            if reaction["type"] == ReactionType.FLOWOUT: 
                indexed_products.append(None)
                indexed_reagents.append(loadedSpecies.index(reaction["in"][0]))

        else: 
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

def excelInit (chemicalSpecies, allParameters, currentTime, refName): 

    chi, delta, ro, Da, div = allParameters[0]
    nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = allParameters[1]
    reactions = allParameters[2]

    #* path directory definition
    directory_name = "../out"

    if not os.path.exists(directory_name):
        try:
            os.makedirs(directory_name)
        except subprocess.CalledProcessError as e:
            print(f"Error in creating the directory: {e}")

    currentData = datetime.now().strftime("%d.%m")
    directory_name = f"../out/out {currentData}"

    if not os.path.exists(directory_name):
        try:
            os.makedirs(directory_name)
        except subprocess.CalledProcessError as e:
            print(f"Error in creating second directory: {e}")

    directory_name = f"../out/out {currentData}/out {currentTime}"

    if not os.path.exists(directory_name):
        try:
            os.makedirs(directory_name)
        except subprocess.CalledProcessError as e:
            print(f"Error in creating second directory: {e}")

    
    name = f"../out/{directory_name}/{refName[1]}.xlsx"
    workbook = xlsxwriter.Workbook (name)

    #* global settings of excel export
    we = workbook.add_worksheet("Environment")
    wk = workbook.add_worksheet("Chemical Species")
    wr = workbook.add_worksheet("Reactions")
    wq = workbook.add_worksheet("Quantity")
    wc = workbook.add_worksheet("Concentration")

    if nFlux > 0:
        wf = workbook.add_worksheet("Flux")

    wq.set_column('A:Z', 15)
    wc.set_column('A:Z', 15)
    we.set_column('A:Z', 18)
    wk.set_column('A:Z', 15)
    wr.set_column('A:Z', 15)
    wr.set_column('B:B', 5)
    
    if nFlux > 0:
        wf.set_column('A:Z', 15)

    even_format = workbook.add_format({'bg_color': '#DDEBF7', 'align': 'center', 'valign': 'vcenter'})
    odd_format = workbook.add_format({'bg_color': '#FFFFFF', 'align': 'center', 'valign': 'vcenter'})
    even_format.set_border(1)
    odd_format.set_border(1)
    
    header_format = workbook.add_format({'bg_color': '#008E3E', 'bold': True, 'align': 'center'})
    header_format.set_border(1)

    #* export parameters
    data = {
        "\u03b4": delta,
        "\u03c1": ro,
        "\u03c7": chi,
        "Da": Da,
        "Duplication Threshold": div,
        "flux": nFlux,
        "gen. to expand": (gen_exp:=gen_exp+1),
        "iterations": nIterates,
        "calving": calving,
        "end time": t_end,
        "max step": max_step,
        "min toll.": toll_min,
        "max toll.": toll_max,
    }

    we.write(0, 0, "Parameter", header_format)
    we.write(0, 1, "Value", header_format)
    
    row = 1
    for key, value in data.items():
        
        cell_format = even_format if row % 2 == 0 else odd_format
        we.write(row, 0, key, cell_format)
        we.write(row, 1, value, cell_format)

        row += 1

    we.set_column('C:XFD', None, None, {'hidden': True})
    we.set_default_row(hide_unused_rows=True)
    
    #* export chemical species
    wk.write(0, 0, "Species", header_format)
    wk.write(0, 1, "Quantity [KG]", header_format)
    wk.write(0, 2, "Coefficient", header_format)

    row = 1
    for i, (species, (quantity, coefficient)) in enumerate(chemicalSpecies.items(), start=1):
        
        cell_format = even_format if row % 2 == 0 else odd_format
        wk.write(row, 0, species, cell_format)
        wk.write(row, 1, quantity, cell_format)
        wk.write(row, 2, coefficient, cell_format)

        row += 1

    wk.set_column('D:XFD', None, None, {'hidden': True})
    wk.set_default_row(hide_unused_rows=True)

    #* export reactions
    _header_format = workbook.add_format({'bg_color': '#FFFF00', 'bold': True, 'align': 'center'})
   
    wr.write(0, 0, "Reagents", header_format)
    wr.write(0, 1, ">", _header_format)
    wr.write(0, 2, "Products", header_format)
    wr.write(0, 3, "Kinetic Coefficient", header_format)
    wr.write(0, 4, "Type", header_format)

    row = 1
    for i, reaction in enumerate(reactions, start=1):
        
        cell_format = even_format if row % 2 == 0 else odd_format
        
        reagents_str = " + ".join(reaction["in"])
        products_str = " + ".join(reaction["out"])

        wr.write(row, 0, reagents_str, cell_format)
        wr.write(row, 1, ">", _header_format)
        wr.write(row, 2, products_str, cell_format)
        wr.write(row, 3, reaction['k'], cell_format)
        wr.write(row, 4, reaction['type'].value, cell_format)
        
        row += 1

    wr.set_column('F:XFD', None, None, {'hidden': True})
    wr.set_default_row(hide_unused_rows=True)

    #* Header writing export
    loadedSpecies = list(chemicalSpecies.keys())
    
    if refName[0] == 1:
        wq.write(0, 0, "Iterations", header_format)
        wc.write(0, 0, "Iterations", header_format)
        wq.write(0,1,"Time", header_format)
        wc.write(0,1,"Time", header_format)

        i = 2
        for species in loadedSpecies: 
            wq.write(0, i, species, header_format)
            wc.write(0, i, species, header_format)
            i+=1
    
    else: 
        wq.write(0, 0, "Generation", header_format)
        wc.write(0, 0, "Generation", header_format)

        i = 1
        for species in loadedSpecies: 
            wq.write(0, i, species, header_format)
            wc.write(0, i, species, header_format)
            i+=1

        wq.write(0, i, "Time", header_format)
        wc.write(0, i, "Time", header_format)

    #wq.set_column(i, 16383, None, {'hidden': True})
    #wc.set_column(f'{i}:{16383}', None, None, {'hidden': True})

    #* writing header about fluxes
    if nFlux > 0:   
        
        ic = 0

        # Expand generation file: [iteration + time + start]
        if refName[0] == 1: 
            wf.write(0, 0, "Iteration", header_format)
            wf.write(0, 1, "Time", header_format)
            ic = 2

            for i, reaction in enumerate(reactions[:nFlux]):
            
                ind = ic
                url = f'internal:Reactions!{ind}:{ind}'
                tx = f"Flux {i+1}"
                
                wf.write_url(0, ic, url, string=tx, cell_format = header_format)

                reagents_str = " + ".join(reaction["in"])
                products_str = " + ".join(reaction["out"])
                
                if reaction["type"] is not ReactionType.FLOWIN and reaction["type"] is not ReactionType.FLOWOUT:
                    wf.write_comment(0, i:=i+2, f'{reagents_str} -> {products_str}')
                
                else: 
                    arrow =" \u2192 "
                    if reaction["type"] is ReactionType.FLOWIN: 
                        wf.write_comment(0, i:=i+2, f'{products_str} {arrow} [CSTR]')
                    else: 
                        wf.write_comment(0, i:=i+2, f'[CSTR] {arrow} {reagents_str}')
                
                ic+=1

        # All generations file: [generation + start + time]
        else: 
            wf.write(0, 0, "Generation", header_format)
            wf.write(0, nFlux+1, "Time", header_format)
            ic = 1

            for i, reaction in enumerate(reactions[:nFlux]):
            
                ind = ic + 1
                url = f'internal:Reactions!{ind}:{ind}'
                tx = f"Flux {i+1}"
                
                wf.write_url(0, ic, url, string=tx, cell_format = header_format)

                reagents_str = " + ".join(reaction["in"])
                products_str = " + ".join(reaction["out"])
                
                if reaction["type"] is not ReactionType.FLOWIN and reaction["type"] is not ReactionType.FLOWOUT:
                    wf.write_comment(0, i:=i+1, f'{reagents_str} -> {products_str}')
                
                else: 
                    arrow =" \u2192 "
                    if reaction["type"] is ReactionType.FLOWIN: 
                        wf.write_comment(0, i:=i+1, f'{products_str} {arrow} [CSTR]')
                    else: 
                        wf.write_comment(0, i:=i+1, f'[CSTR] {arrow} {reagents_str}')
                
                ic+=1

    if nFlux == 0:
        return workbook, wc, wq
    
    if nFlux > 0:
        return workbook, wc, wq, wf

def excelExport1 (matrixSimulation, timeSimulation, chemicalSpecies, allParameters, refName): 

    #* path directory definition
    directory_name = "../out"

    if not os.path.exists(directory_name):
        try:
            os.makedirs(directory_name)
        except subprocess.CalledProcessError as e:
            print(f"Error in creating the directory: {e}")

    currentData = datetime.now().strftime("%d.%m")
    directory_name = f"../out/out {currentData}"

    if not os.path.exists(directory_name):
        try:
            os.makedirs(directory_name)
        except subprocess.CalledProcessError as e:
            print(f"Error in creating second directory: {e}")

    currentTime = datetime.now().strftime("%H.%M.%S")
    name = f"../out/{directory_name}/{refName[1]} {currentTime}.xlsx"
    workbook = xlsxwriter.Workbook (name)

    #* global settings of excel export
    we = workbook.add_worksheet("Environment")
    wk = workbook.add_worksheet("Chemical Species")
    wr = workbook.add_worksheet("Reactions")
    wq = workbook.add_worksheet("Quantity")
    wc = workbook.add_worksheet("Concentration")

    wq.set_column('A:Z', 15)
    wc.set_column('A:Z', 15)
    we.set_column('A:Z', 18)
    wk.set_column('A:Z', 15)
    wr.set_column('A:Z', 15)
    wr.set_column('B:B', 5)

    even_format = workbook.add_format({'bg_color': '#DDEBF7', 'align': 'center', 'valign': 'vcenter'})
    odd_format = workbook.add_format({'bg_color': '#FFFFFF', 'align': 'center', 'valign': 'vcenter'})
    even_format.set_border(1)
    odd_format.set_border(1)
    
    header_format = workbook.add_format({'bg_color': '#008E3E', 'bold': True, 'align': 'center'})
    header_format.set_border(1)

    #* export parameters
    chi, delta, ro, Da, div = allParameters[0]
    nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = allParameters[1]
    reactions = allParameters[2]

    data = {
        "\u03b4": delta,
        "\u03c1": ro,
        "\u03c7": chi,
        "Da": Da,
        "Duplication Threshold": div,
        "flux": nFlux,
        "gen. to expand": (gen_exp:=gen_exp+1),
        "iterations": nIterates,
        "calving": calving,
        "end time": t_end,
        "max step": max_step,
        "min toll.": toll_min,
        "max toll.": toll_max,
    }

    we.write(0, 0, "Parameter", header_format)
    we.write(0, 1, "Value", header_format)
    
    row = 1
    for key, value in data.items():
        
        cell_format = even_format if row % 2 == 0 else odd_format
        we.write(row, 0, key, cell_format)
        we.write(row, 1, value, cell_format)

        row += 1

    we.set_column('C:XFD', None, None, {'hidden': True})
    we.set_default_row(hide_unused_rows=True)
    
    #* export chemical species
    wk.write(0, 0, "Species", header_format)
    wk.write(0, 1, "Quantity [KG]", header_format)
    wk.write(0, 2, "Coefficient", header_format)

    row = 1
    for i, (species, (quantity, coefficient)) in enumerate(chemicalSpecies.items(), start=1):
        
        cell_format = even_format if row % 2 == 0 else odd_format
        wk.write(row, 0, species, cell_format)
        wk.write(row, 1, quantity, cell_format)
        wk.write(row, 2, coefficient, cell_format)

        row += 1

    wk.set_column('D:XFD', None, None, {'hidden': True})
    wk.set_default_row(hide_unused_rows=True)

    #* export reactions
    _header_format = workbook.add_format({'bg_color': '#FFFF00', 'bold': True, 'align': 'center'})
   
    wr.write(0, 0, "Reagents", header_format)
    wr.write(0, 1, ">", _header_format)
    wr.write(0, 2, "Products", header_format)
    wr.write(0, 3, "Kinetic Coefficient", header_format)
    wr.write(0, 4, "Type", header_format)

    row = 1
    for i, reaction in enumerate(reactions, start=1):
        
        cell_format = even_format if row % 2 == 0 else odd_format
        
        reagents_str = " + ".join(reaction["in"])
        products_str = " + ".join(reaction["out"])

        wr.write(row, 0, reagents_str, cell_format)
        wr.write(row, 1, ">", _header_format)
        wr.write(row, 2, products_str, cell_format)
        wr.write(row, 3, reaction['k'], cell_format)
        wr.write(row, 4, reaction['type'].value, cell_format)
        
        row += 1

    wr.set_column('F:XFD', None, None, {'hidden': True})
    wr.set_default_row(hide_unused_rows=True)

    #* export quantity and concentration info
    
    # Header writing
    if refName[0] == 1:
        wq.write(0, 0, "Iterations", header_format)
        wc.write(0, 0, "Iterations", header_format)
    else: 
        wq.write(0, 0, "Generation", header_format)
        wc.write(0, 0, "Generation", header_format)
    
    # Index writing
    for i in range(1, len(matrixSimulation) + 1):
        wq.write(i,0,i)
        wc.write(i,0,i)

    wq.write(0,1,"Time", header_format)
    wc.write(0,1,"Time", header_format)

    loadedSpecies = list(chemicalSpecies.keys())
    
    i = 2
    for species in loadedSpecies: 
        wq.write(0, i, species, header_format)
        wc.write(0, i, species, header_format)
        i+=1

    #wq.set_column(i, 16383, None, {'hidden': True})
    #wc.set_column(f'{i}:{16383}', None, None, {'hidden': True})

    # Timing export
    row=1 
    column = 1
    for value in np.array(timeSimulation):
        wq.write(row,column,value)
        wc.write(row,column,value)
        row+=1

    # Data export [quantities and concentrations]
    row=1
    for matLine in matrixSimulation:
        column=2
        for value in matLine:
            wq.write(row, column, value)
            wc.write(row, column, getConcentration (value, matLine[0], allParameters[0][2], allParameters[0][1]))
            column+=1
        row+=1

    wc.set_default_row(hide_unused_rows=True)
    wq.set_default_row(hide_unused_rows=True)

    workbook.close()

def excelExport (matrixSimulation, timeSimulation, chemicalSpecies, allParameters, currentTime, refName): 

    nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = allParameters[1]

    if nFlux == 0:
        workbook, wc, wq = excelInit(chemicalSpecies, allParameters, currentTime, refName)

    if nFlux > 0:
        workbook, wc, wq, wf = excelInit(chemicalSpecies, allParameters, currentTime, refName)

    #* Index writing
    if nFlux == 0:
        for i in range(1, len(matrixSimulation) + 1):
            wq.write(i,0,i)
            wc.write(i,0,i)

    if nFlux > 0:
        for i in range(1, len(matrixSimulation) + 1):
            wq.write(i,0,i)
            wc.write(i,0,i)
            wf.write(i,0,i)

    #* Timing export
    if refName[0] == 1: 
        row=1 
        column = 1

        if nFlux == 0:
            for value in np.array(timeSimulation):
                wq.write(row,column,value)
                wc.write(row,column,value)
                row+=1
    
        if nFlux > 0:
            for value in np.array(timeSimulation):
                wq.write(row,column,value)
                wc.write(row,column,value)
                wf.write(row,column,value)
                row+=1
    
    else: 
        row=1 
        column = len(chemicalSpecies)
        
        if nFlux == 0:
            for value in np.array(timeSimulation):
                wq.write(row,column,value)
                wc.write(row,column,value)
                row+=1


        if nFlux > 0:
            for value in np.array(timeSimulation):
                wq.write(row,column,value)
                wc.write(row,column,value)
                row+=1

            row=1 
            column = nFlux+1
            for value in np.array(timeSimulation):
                wf.write(row,column,value)
                row+=1 
        
        
    #* Data export [quantities and concentrations] + [any export fluxes]

    chemicalVariation = [row[:len(chemicalSpecies)] for row in matrixSimulation]
    if nFlux > 0: 
        fluxes = [row[-nFlux:] for row in matrixSimulation]

    if refName[0] == 1: 
        
        if nFlux == 0:
            row=1
            for matLine in chemicalVariation:
                column=2
                for value in matLine:
                    wq.write(row, column, value)
                    wc.write(row, column, getConcentration (value, matLine[0], allParameters[0][2], allParameters[0][1]))
                    column+=1
                row+=1
        
        if nFlux > 0: 
            row=1
            for matLine in chemicalVariation:
                column=2
                for value in matLine:
                    wq.write(row, column, value)
                    wc.write(row, column, getConcentration (value, matLine[0], allParameters[0][2], allParameters[0][1]))
                    column+=1
                row+=1

            row=1
            for matLine in fluxes:
                column=2
                for value in matLine:
                    wf.write(row, column, value)
                    column+=1
                row+=1

    else: 
        if nFlux == 0:
            row=1
            for matLine in chemicalVariation:
                column=1
                for value in matLine:
                    wq.write(row, column, value)
                    wc.write(row, column, getConcentration (value, matLine[0], allParameters[0][2], allParameters[0][1]))
                    column+=1
                row+=1
        
        if nFlux > 0:
            row=1
            for matLine in chemicalVariation:
                column=1
                for value in matLine:
                    wq.write(row, column, value)
                    wc.write(row, column, getConcentration (value, matLine[0], allParameters[0][2], allParameters[0][1]))
                    column+=1
                row+=1

            row=1
            for matLine in fluxes:
                column=1
                for value in matLine:
                    wf.write(row, column, value)
                    column+=1
                row+=1
            
    wc.set_default_row(hide_unused_rows=True)
    wq.set_default_row(hide_unused_rows=True)
    
    if nFlux > 0: 
        wf.set_default_row(hide_unused_rows=True)

    workbook.close()

def printFinalInfo (currentTime, parameters, environment, chemicalSpecies, reactions, matrixSimulation): 

    print ("\n->END SIMULATION<-")
    printInfo(parameters, environment, chemicalSpecies, reactions)

    print ("\n->About last generation<-")
    lastGen = matrixSimulation[-1]

    i=0
    for species in chemicalSpecies.items():
        print (f"{i+1}.", species[0], f"[{lastGen[i]} Kg]\n")
        i+=1
    
    print ("\n->About flux in last generation<-")
    nFlux = environment[5]
    
    for i in range (len(chemicalSpecies), len(lastGen)):
        print(f"{i+1-len(chemicalSpecies)}. Flux", f"[{lastGen[i]}]\n")
        i+=1

    print (f"All simulation data has been successfully exported to the directory /out/sim {currentTime}.")