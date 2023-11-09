from enum import Enum

from errorsCheck import checkProtoSim

class ReactionType (Enum): 
    
    CONDENSATION_21 = "Condensation: II>O"
    CONDENSATION_22 = "Condensation: II>OO"

    CLEAVAGE_12 = "Cleavage: I>OO"
    CLEAVAGE_23 = "Cleavage: II>OOO"
    
    DIFFUSION = "Diffusion"

    ND = "Undefined Type"

def identifyType (reaction, verbose): 

    reactionsParts = reaction.split(">")
   
    reactants = reactionsParts[0].strip()
    products = reactionsParts[1].strip()
    
    nReactans = len (reactants.split ("+"))
    nProducts = len (products.split ("+"))

    if verbose: 
        print(f"Identifying Reaction Type: {nReactans};{nProducts}")

    if nReactans == 1: 
        if nProducts == 1: 
            return ReactionType.DIFFUSION
        if nProducts == 2: 
            return ReactionType.CLEAVAGE_12

    if nReactans == 2: 
        if nProducts == 1: 
            return ReactionType.CONDENSATION_21
        if nProducts == 2: 
            return ReactionType.CONDENSATION_22
        if nProducts == 3: 
            return ReactionType.CLEAVAGE_23

    checkProtoSim(2, reaction)