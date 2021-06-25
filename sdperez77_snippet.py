    # new code for trait probabilities
  for person in people:
        ngenes = [person]['ngenes']
        if person in have_trait:
    # Fourth scenario - person in have_trait
            probdict[person]['prob4'] = PROBS['trait'][ngenes][True]
    # Fifth scenario is everyone not in the set have_trait doesn't have the trait
        else:
            probdict[person]['prob5'] = PROBS['trait'][ngenes][False]
        
    # Multiply joint probabilities for all scenarios for all people to determine joint probability
    # modified to remove cuplicate for loops and use *= operator
    for person in probdict:
        if probdict[person]['prob1'] != 0:
            jp *= probdict[person]['prob1']
        if probdict[person]['prob2'] != 0:
            jp *= probdict[person]['prob2']
        if probdict[person]['prob3'] != 0:
            jp *= probdict[person]['prob3']
        if probdict[person]['prob4'] != 0:
            jp *= probdict[person]['prob4']
        if probdict[person]['prob5'] != 0:
            jp *= probdict[person]['prob5']
