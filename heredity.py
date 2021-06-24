import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    # Create empty dictionary to store probabilities for each person for each scenario, setting all probabilities initially equal to 0.
    probdict = {
        person: {
            'prob1': 0,
            'prob2': 0, 
            'prob3': 0,
            'prob4': 0,
            'prob5': 0 
        }
        for person in people
    }

    # Create variable joint probability, jp, initially set to 1 (since we will be multiplying iteratively)
    jp = 1
    # First scenario is probabilities for everyone in set one_gene has one copy of gene
    for person in one_gene:
        # If the person has no mother listed, then probability will be just the distribution in the population in general
        if people[person]['mother'] == None and people[person]['father'] == None:
            probdict[person]['prob1'] = PROBS['gene'][1]
        # If the person does have a mother listed, then probability will be combination of receiving the gene from mother and father
        elif people[person]['mother'] != None and people[person]['father'] != None:
            if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                # Probability of receiving gene from either parent is solely based on mutation given neither parent has the gene
                probparent = 0.01
                # Probability of not receiving gene from either parent is solely based on mutation given neither parent has the gene
                probnotparent = 0.99
                probdict[person]['prob1'] = 2 * probparent * probnotparent
            if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] in one_gene:
                # Probability of receiving gene from mother is solely based on mutation given mother does not have the gene
                probM = 0.01
                probnotM = 0.99
                # Probability of receiving gene from father is 50% chance of not getting the gene (unless mutated) + 50% of getting the gene (unless mutated)
                probF = (0.50 * 0.01 + 0.50 * 0.99)
                probnotF = (0.50 * 0.99 + 0.50 * 0.01)
                probdict[person]['prob1'] = probM * probnotF + probF * probnotM
            if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] in two_genes:
                probM = 0.01
                probnotM = 0.99
                probF = 0.99
                probnotF = 0.01
                probdict[person]['prob1'] = probM * probnotF + probF * probnotM
            if people[person]['mother'] in one_gene and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                probM = (0.50 * 0.01 + 0.50 * 0.99)
                probnotM = (0.50 * 0.99 + 0.50 * 0.01)
                probF = 0.01
                probnotF = 0.99
                probdict[person]['prob1'] = probM * probnotF + probF * probnotM
            if people[person]['mother'] in one_gene and people[person]['father'] in one_gene:
                probparent = (0.50 * 0.01 + 0.50 * 0.99)
                probnotparent = (0.50 * 0.99 + 0.50 * 0.01)
                probdict[person]['prob1'] = 2 * probparent * probnotparent 
            if people[person]['mother'] in one_gene and people[person]['father'] in two_genes:
                probM = (0.50 * 0.01 + 0.50 * 0.99)
                probnotM = (0.50 * 0.99 + 0.50 * 0.01)
                probF = 0.99
                probnotF = 0.01
                probdict[person]['prob1'] = probM * probnotF + probF * probnotM
            if people[person]['mother'] in two_genes and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                probM = 0.99
                probnotM = 0.01
                probF = 0.01
                probnotF = 0.99
                probdict[person]['prob1'] = probM * probnotF + probF * probnotM
            if people[person]['mother'] in two_genes and people[person]['father'] in one_gene:
                probM = 0.99
                probnotM = 0.01
                probF = (0.50 * 0.01 + 0.50 * 0.99)
                probnotF = (0.50 * 0.99 + 0.50 * 0.01)
                probdict[person]['prob1'] = probM * probnotF + probF * probnotM
            if people[person]['mother'] in two_genes and people[person]['father'] in two_genes:
                probparent = 0.99
                probnotparent = 0.01
                probdict[person]['prob1'] = 2 * probparent * probnotparent
    # Second scenario is probabilities for everyone in set two_genes have two copies of gene
    for person in two_genes:
        if people[person]['mother'] == None and people[person]['father'] == None:
            probdict[person]['prob2'] = PROBS['gene'][2]
        elif people[person]['mother'] != None and people[person]['father'] != None:
            if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                probM = 0.01
                probF = 0.01
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] in one_gene:
                probM = 0.01
                probF = (0.50 * 0.01 + 0.50 * 0.99)
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] in two_genes:
                probM = 0.01
                probF = 0.99
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] in one_gene and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                probM = (0.50 * 0.01 + 0.50 * 0.99)
                probF = 0.01
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] in one_gene and people[person]['father'] in one_gene:
                probM = (0.50 * 0.01 + 0.50 * 0.99)
                probF = (0.50 * 0.01 + 0.50 * 0.99)
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] in one_gene and people[person]['father'] in two_genes:
                probM = (0.50 * 0.01 + 0.50 * 0.99)
                probF = 0.99
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] in two_genes and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                probM = 0.99
                probF = 0.01
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] in two_genes and people[person]['father'] in one_gene:
                probM = 0.99
                probF = (0.50 * 0.01 + 0.50 * 0.99)
                probdict[person]['prob2'] = probM * probF
            if people[person]['mother'] in two_genes and people[person]['father'] in two_genes:
                probM = 0.99
                probF = 0.99
                probdict[person]['prob2'] = probM * probF
    # Third scenario is probabilities for everyone who is not in either one gene or two gene sets does not have a copy of the gene:
    for person in people:
        if person not in one_gene and person not in two_genes:
            if people[person]['mother'] == None and people[person]['father'] == None:
                probdict[person]['prob3'] = PROBS['gene'][0]
            elif people[person]['mother'] != None and people[person]['father'] != None:
                if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                    probM = 0.99
                    probF = 0.99
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] in one_gene:
                    probM = 0.99
                    probF = (0.50 * 0.01 + 0.50 * 0.99)
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] not in one_gene and people[person]['mother'] not in two_genes and people[person]['father'] in two_genes:
                    probM = 0.99
                    probF = 0.01
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] in one_gene and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                    probM = (0.50 * 0.01 + 0.50 * 0.99)
                    probF = 0.99
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] in one_gene and people[person]['father'] in one_gene:
                    probM = (0.50 * 0.01 + 0.50 * 0.99)
                    probF = (0.50 * 0.01 + 0.50 * 0.99)
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] in one_gene and people[person]['father'] in two_genes:
                    probM = (0.50 * 0.01 + 0.50 * 0.99)
                    probF = 0.01
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] in two_genes and people[person]['father'] not in one_gene and people[person]['father'] not in two_genes:
                    probM = 0.01
                    probF = 0.99
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] in two_genes and people[person]['father'] in one_gene:
                    probM = 0.01
                    probF = (0.50 * 0.01 + 0.50 * 0.99)
                    probdict[person]['prob3'] = probM * probF
                if people[person]['mother'] in two_genes and people[person]['father'] in two_genes:
                    probM = 0.01
                    probF = 0.01
                    probdict[person]['prob3'] = probM * probF
    # Fourth scenario is everyone in the set have_trait has the trait
    for person in have_trait:
        # If the person has one gene then probability of having the trait is given by population distribution
        if person in one_gene:
            probdict[person]['prob4'] = PROBS['trait'][1][True]
        if person in two_genes:
            probdict[person]['prob4'] = PROBS['trait'][2][True]
        if person not in one_gene and person not in two_genes:
            probdict[person]['prob4'] = PROBS['trait'][0][True]
    # Fifth scenario is everyone not in the set have_trait doesn't have the trait
    for person in people:
        if person not in have_trait:
            # If the person has one gene then probability of not having the trait is given by population distribution
            if person in one_gene:
                probdict[person]['prob5'] = PROBS['trait'][1][False]
            if person in two_genes:
                probdict[person]['prob5'] = PROBS['trait'][2][False]
            if person not in one_gene and person not in two_genes:
                probdict[person]['prob5'] = PROBS['trait'][0][False]
    # Multiply joint probabilities for all scenarios for all people to determine joint probability
    for person in probdict:
        if probdict[person]['prob1'] != 0:
            jp = jp * probdict[person]['prob1']
    for person in probdict:
        if probdict[person]['prob2'] != 0:
            jp = jp * probdict[person]['prob2']
    for person in probdict:
        if probdict[person]['prob3'] != 0:
            jp = jp * probdict[person]['prob3']
    for person in probdict:
        if probdict[person]['prob4'] != 0:
            jp = jp * probdict[person]['prob4']
    for person in probdict:
        if probdict[person]['prob5'] != 0:
            jp = jp * probdict[person]['prob5']
    # Return value is joint probability of the scenario
    return jp


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    people = load_data(sys.argv[1])
    
    # Loop over all people in dictionary
    for person in people:
        # If the person is in the one_gene set, then add joint probability p to the joint probability distribution for having one gene:
        if person in one_gene:
            probabilities[person]['gene'][1] = probabilities[person]['gene'][1] + p

    for person in people:
        if person in two_genes:
            probabilities[person]['gene'][2] = probabilities[person]['gene'][2] + p

    for person in people:
        if person in have_trait:
            probabilities[person]['trait'][True] = probabilities[person]['trait'][True] + p

    for person in people:
        if person not in one_gene and person not in two_genes:
            probabilities[person]['gene'][0] = probabilities[person]['gene'][0] + p

    for person in people:
        if person not in have_trait:
            probabilities[person]['trait'][False] = probabilities[person]['trait'][False] + p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    people = load_data(sys.argv[1])
    
    # Create empty dictionary for keeping track of the sum of the values of genes for each person
    genesum = dict()

    # Create empty dictionary for keeping track of the sume of the values of traits for each person
    traitsum = dict()

    # Loop over all people in the dictionary
    for person in people:
        # Calculate sum of genes
        genesum[person] = probabilities[person]['gene'][2] + probabilities[person]['gene'][1] + probabilities[person]['gene'][0]
        # Update probability distribution for each gene so that it is normalized to 1
        probabilities[person]['gene'][2] = probabilities[person]['gene'][2] / genesum[person]
        probabilities[person]['gene'][1] = probabilities[person]['gene'][1] / genesum[person]
        probabilities[person]['gene'][0] = probabilities[person]['gene'][0] / genesum[person]
    for person in people:
        # Calculate sum of traits
        traitsum[person] = probabilities[person]['trait'][True] + probabilities[person]['trait'][False]
        # Update probability distribution for each trait so that it is normalized to 1
        probabilities[person]['trait'][True] = probabilities[person]['trait'][True] / traitsum[person]
        probabilities[person]['trait'][False] = probabilities[person]['trait'][False] / traitsum[person]


if __name__ == "__main__":
    main()
