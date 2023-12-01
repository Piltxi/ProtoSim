import numpy as np
import time
from tqdm import tqdm

from chemicalio import excelExport, excelInit, map_species_to_indices, getConcentration
from reactions import ReactionType
from errorsCheck import checkProtoSim

"""
parameters = allParameters[0]
#parameters = [chi, delta, ro, k, Da, As, div]
chi, delta, ro, k, Da, As, div = parameters

environment = allParameters [1]
# environment = [nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp]; 
nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = environment
"""

def scalar_multiply(vector, scalar):
    result = []
    
    for i in range (len(vector)):
        result += [vector[i]*scalar]

    return result

def add_vectors(vector1, vector2):
    result = []
    
    for i in range(len(vector1)):
        result += [vector1[i] + vector2[i]]
    
    return result

def divisionTest(time, protoAct, parameters):

    if protoAct[0] > parameters[-1]:
        return False
    
    else:
        return True

def tolleranceTest(protoAct, protoNext, s_min, s_max, nFlussi, dt):

    if not protoNext:
        return False
        
    for i in range(len(protoAct)-nFlussi):
       
        if protoNext[i] < 0:
            return True
        
        if protoNext[i]!=0 and protoAct[i] != 0:
            var = protoNext[i] / protoAct[i]
            if var > s_max or var < s_min:
                return True
    
    return False

def callOdeSolver (ode_function, time, protoAct, parameters, mapReactions, deltaT, nFlux):

    for i in range(len(protoAct[1]) - nFlux):
        if protoAct[1][i]<0:
            checkProtoSim (4, [protoAct[1][i], i])
    
    var = ode_function (time, protoAct, [parameters, mapReactions])

    return var

def ode_function (time, protoAct, parameters): 

    protoX = protoAct[1][:]
    coefficients = protoAct[0] [:]
    
    Dx=scalar_multiply(protoX, 0)


    """
        Primi 
    """

    for i in range(len(protoX)):
        if coefficients[i]!=0:
            Dx[0] += protoX[i] * coefficients[i]

    reactions = parameters[1]
    
    # Reaction Variation Rules
    for i in range(len(reactions)):
        
        match reactions[i]["type"]:

            case ReactionType.FLOWIN: 
                
                protoX[reactions[i]["out"][0]] =+ reactions[i]["k"]

            case ReactionType.FLOWOUT: 

                term = reactions[i]["k"] * protoX[reactions[i]["in"][0]]
                protoX[reactions[i]["in"][0]] =- term

            case ReactionType.CONDENSATION_21: 
                
                term = reactions[i]["k"] * protoX[reactions[i]["in"][0]] * protoX [reactions[i]["in"][1]] / (parameters[0][0] * pow (protoX[0], 1.5))
                Dx[reactions[i]["in"][0]] -= term
                Dx[reactions[i]["in"][1]] -= term
                Dx[reactions[i]["out"][0]] += term
           
            case ReactionType.CONDENSATION_22: 
    
                term = reactions[i]["k"] * protoX[reactions[i]["in"][0]] * protoX[reactions[i]["in"][1]] / (parameters[0][0] * pow (protoX[0], 1.5))
                Dx[reactions[i]["in"][0]] -= term
                Dx[reactions[i]["in"][1]] -= term
                Dx[reactions[i]["out"][0]] += term
                Dx[reactions[i]["out"][1]] += term
            
            case ReactionType.CLEAVAGE_12: 
                
                term = reactions[i]["k"] * protoX[reactions[i]["in"][0]]
                Dx[reactions[i]["in"][0]] -= term
                Dx[reactions[i]["out"][0]] += term
                Dx[reactions[i]["out"][1]] += term
            
            case ReactionType.CLEAVAGE_23: 
                
                term = reactions[i]["k"] * protoX[reactions[i]["in"][0]] * protoX[reactions[i]["in"][1]] / (parameters[0][0] * pow (protoX[0], 1.5))
                Dx[reactions[i]["in"][0]] -= term
                Dx[reactions[i]["in"][1]] -= term
                Dx[reactions[i]["out"][0]] += term
                Dx[reactions[i]["out"][1]] += term
                Dx[reactions[i]["out"][2]] += term

            case ReactionType.DIFFUSION:

                Dx [reactions[i]["out"][0]] += ((parameters[0][4] * ( protoX[0] / (parameters [0][2] * parameters[0][1]) ) * parameters[0][3]) * (reactions[i]["in"][0] - (protoX[reactions[i]["out"][0]] / (parameters[0][0] * pow (protoX[0], 1.5))))) / parameters[0][1]

            case _:
                checkProtoSim (5, reactions[i])

    return Dx

def solver (ode_function, interval, protoGen, mapReactions, parameters, environment, divisionTest, maxStep, tollerance, nFlux, coefficient):

    deltaT = min (maxStep/10., interval[1]/10.)
    
    # Start of current simulation
    t = interval[0]
    
    protoAct = protoGen[1][:] 

    # Resolution of negative quantities
    for i in range (len(protoAct)-nFlux):
            if protoAct[i]<0: 
                    protoAct[i] = 0 
 
    tempi = []
    y = []    
    
    tempi += [t]
    y += [protoAct[:]]
    
    seconds = 0.01

    while divisionTest (t, protoAct, parameters):
        
        # if t > seconds: 
        #     print (t, protoAct)
        
        var = callOdeSolver (ode_function, t, [coefficient, protoAct], parameters, mapReactions, deltaT, nFlux)
        
        protoNext = add_vectors (protoAct, scalar_multiply(var, deltaT))

        if not tolleranceTest (protoAct, protoNext, tollerance[0], tollerance[1], nFlux, deltaT):
            deltaT *= 1.2
            if deltaT > maxStep:
                deltaT = maxStep
            
        while tolleranceTest (protoAct, protoNext, tollerance[0], tollerance[1], nFlux, deltaT):
            deltaT /= 2
            protoNext = add_vectors (protoAct, scalar_multiply(var, deltaT))

        if t+deltaT > interval [1]:
            deltaT = interval[1]-t

        t += deltaT
        tempi += [t]
        y += [protoNext]

        protoAct = protoNext [:]

        if t >= interval [1]: 
            print ("End; ", t, interval[1])
            break
    
    return (tempi, y)

def simulation (verbose, environment, parameters, chemicalSpecies, reactions): 

    nIterates, t_end, max_step, toll_min, toll_max, nFlux, gen_exp, calving = environment

    # Preparing to directly export data to the csv file.
    workbook, wc, wq = excelInit(chemicalSpecies, [parameters, environment, reactions])

    protoGen = np.array([chemicalSpecies[quantity][0] for quantity in chemicalSpecies])
    protoInit = np.array([chemicalSpecies[quantity][0] for quantity in chemicalSpecies])
    loadedSpecies = list(chemicalSpecies.keys())
    mapReactions = map_species_to_indices(reactions, loadedSpecies)
    coefficient = np.array ([chemicalSpecies[coefficient][1] for coefficient in chemicalSpecies])

    time_ = []
    mat = []
    
    t_start = 0.

    if not verbose:
        progress_bar = tqdm(total=nIterates, desc="Simulating", unit="generation", position=0, dynamic_ncols=True)

    try: 
        for i in range(nIterates):

            if verbose: 
                print ("Start generation n.", i+1)

            # num_sol = solve_ivp(ode_fn, [t_begin, t_end], [x_init], method=method, dense_output=True)
            startTime = time.time()
            (solverTime, y_sol) = solver (ode_function, [t_start, t_end], [protoInit, protoGen], mapReactions, parameters, environment, divisionTest, max_step, [toll_min, toll_max], nFlux, coefficient)
            endTime = time.time()

            #print ("Duplication Time: ", solverTime[-1])
            
            if i == gen_exp: 
                excelExport(y_sol, solverTime, chemicalSpecies, [parameters, environment, reactions], [1, "expand"])

            executionTime = endTime - startTime
            if verbose: 
                if executionTime > 60:
                    minutes=int(executionTime/60)
                    seconds=round(executionTime%60)
                    print(f"Duplication Time: ", solverTime[-1], f" | Time spent {minutes}:{seconds} minutes")
                else:  
                    print(f"Duplication Time: ", solverTime[-1], f" | time spent {round(executionTime)} seconds") 
            
            # About last generation 
            protoGen = np.copy (y_sol[-1])
        
            if verbose: 
                print ("End generation n.", i+1, "\t", protoGen, "\n")
        
            time_ += [solverTime[-1]]
            mat += [np.copy(protoGen)]

            wq.write(i+1, 0, i+1)
            for j, (species, _) in enumerate(chemicalSpecies.items(), start=1):
                wq.write(i+1, j, protoGen[j-1])
            wq.write(i + 1, len(chemicalSpecies) + 1, solverTime[-1])

            wc.write(i+1, 0, i+1)
            for j, (species, _) in enumerate(chemicalSpecies.items(), start=1):
                wc.write(i+1, j, getConcentration (protoGen[j-1], protoGen[0], parameters[2], parameters[1]))
            wc.write(i + 1, len(chemicalSpecies) + 1, solverTime[-1])

            # Duplication
            protoGen[0] = protoGen[0]/2.
            
            # Without Flux Analysis
            for i in range(1,len(protoGen)-nFlux):
                protoGen[i]=protoGen[i]*calving
            
            # With Flux Analysis
            for i in range(nFlux):
                protoGen[-1-i] = 0 
            
            if not verbose:
                progress_bar.update(1)
    
    except KeyboardInterrupt:
        
        if not verbose:
            progress_bar.close()
        
        workbook.close()
        print ("\nimproper shutdown - current simulation exported successfully")
        quit()

    if not verbose:
        progress_bar.close()

    workbook.close()

    return (time_, mat)
