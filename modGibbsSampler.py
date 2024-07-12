import sys
import random
import numpy


def gibbs_sampler_actual(dna, k: int, t: int, n: int) -> list[str]:
    """Implements the GibbsSampling algorithm for motif finding."""

    # Arbitrarily large function counter
    function_counter = 10
    
    # Call 'get_random_kmers()' by passing in dna and k.  Will return a list of k-length strings.
    random_kmer_list = get_random_kmers(dna, k)

    
    # Benchmark lowest_scoring_kmers. Uses a copy so that the random_kmers_list doesnt get modified in the code. 
    best_kmers = random_kmer_list.copy()

    # Loop the folling actions 1000 times.
    while function_counter > 0:
        # for tracking number of loops. 
        function_counter -= 1

        

        # Call isolate_string_and_kmer_list.  The function will randomly pull one of the dna strings out
        # and will do the same for the kmer list.
        isolated_string, modified_kmer_list, random_position = isolate_and_modify(dna, best_kmers, t)

        # Form a profile matrix based on the modified k-mer list, i.e., the list where
        # a kmer was removed at random.
        profile = form_profile(modified_kmer_list, k)

        #Use the profile to determine probabilities of all kmers in isolated string.
        # Then roll a biased die weighted with those probabilities to select a kmer.
        die_roll_kmer = find_most_probable_kmer(isolated_string, k, profile)
        #print('die roll', die_roll_kmer, 'pos', random_position)

        # Restore the die_roll_kmer to the kmer_list in the correct position.         
        modified_kmer_list.insert(random_position, die_roll_kmer)


        # Score the modified list against the benchmark list:
        score_modified = score(modified_kmer_list, k, t)
        score_best = score(best_kmers, k, t)

        #print(best_kmers, score_best)
        #print(modified_kmer_list, score_modified)

        #print('mod score', score_modified, 'score best', score_best)

        # Compare the score of the algorithmlically modified list to the previous best kmer list. 
        if score_modified < score_best:
            # Whichever has the lowest score is the benchmark list.
            best_kmers = modified_kmer_list

    #print(best_kmers, score_best)
    return best_kmers, score_best


def gibbs_sampler(dna: list[str], k: int, t: int, n: int) -> list[str]:

    """Run the Gibbs sampler multiple times to avoid local minima."""
    best_score, best_motifs = float('inf'), None

    for _ in range(n):
        motifs, score = gibbs_sampler_actual(dna, k, t, n)
        
        if score < best_score and score > 0:
            best_score = score
            best_motifs = motifs

            #print('score benchmark', best_score)
    return best_motifs, best_score



def get_random_kmers(dna: list[str], k: int) -> list[str]:
    """ This function randomly selects kmers in each string of dna"""
    
    # length of a given dna string
    n = len(dna[0])
    kmer_count = n - k - 1
    #print('kmer count', kmer_count)

    random_kmers = []
    
        # Iterate strings of dna
    for i in range(len(dna)):
        # Generate random number between 1 and kmer_count
        random_kmer_int = random.randint(0, kmer_count)

        # Find number of kmers in a given dna string
        # variable for random kmer in the i-th string
        random_kmer = dna[i][random_kmer_int:random_kmer_int + k]
        #print(random_kmer)
            
        # add the variable kmers per list in the i loop, not the j loop. 
        random_kmers.append(random_kmer)

    return random_kmers

def isolate_and_modify(dna, motifs, t):
    """ This function randomly chooses one of the dna strings and isolates it.  It also
    goes into the motifs list and removes the kmer associated with the removed dna string.
    For example, if we isolate the 3rd dna string, we will also remove the 3rd kmer."""

    # This prevents motifs from being modified.
    copy_of_motifs = motifs.copy()

    # Pick a random position between 0 and t.. Must substract one. 
    random_position = random.randint(0, t-1)

    # Get the string associated with a random position in the dna list.  Isolate it.
    isolated_string = dna[random_position]

    # Pop (remove) the kmer from the motifs list with the same random position.
    copy_of_motifs.pop(random_position)

    return isolated_string, copy_of_motifs, random_position

    


def form_profile(motifs, k):
    """Generate the profile matrix for a given list.  Find probability of each bases in
    each column in the kmer list."""
    """Note: this was copy-pasted from 2.6 assigment. It makes use of pseudocounts."""
    
    profile = []
    t = len(motifs)
    # Initialize the profile with zeros
    for _ in range(k):
        profile.append({'A': 1.0, 'C': 1.0, 'G': 1.0, 'T': 1.0})
    # Count each base occurrence in each column
    for i in range(k):  # iterate rows
        for motif in motifs:    # Iterate letters/columns
            
            # profile[i] is the row
            # Motif[i] is the key in each row.
            # By convention of retrieving dictionary value via dictionary[key] = value, we do:
            
            profile[i][motif[i]] += 1
    
    # Divide each iposition in the matrix by t
    for row in profile:
        for key in row:
            row[key] /= (t + 4)


    return profile


def find_most_probable_kmer(string, k, profile):
    """ Based on a matrix of base probabilities and an input dna string, this function will
    go over every kmer in the dna string and find the one with the closest matching probability
    to the matrix"""

    # List of probabilities of each kmer, divided by the sum of probabilities.
    prob_distribution = []
    
    # Baseline sum of probabilities for the dna string.
    prob_sum = 0
    kmer_list = []

    n = len(string)
    # iterate over the dna string and pull out kmers.
    for i in range(n - k + 1):
        kmer = string[i:i + k] # each kmer.
        kmer_list.append(kmer)
        prob = 1    # Default probability will be 1 so we can multiply.
        # for each position in kmer, and each character in kmer (each base),  multiply the probability
        # against the corresponding probability in our matrix profile for that position.
        # store the new probability value in 'prob', our new benchmark.

        #print('kmer', kmer)    # Test print
        
        for j, char in enumerate(kmer):
            
            prob *= profile[j][char]
            # Iteratively sort our highest probability and store in max_prob.

        # once a kmer is assigned a probability, chuck the value into a sum equation
        prob_sum += prob

        # Add each k-mer's probability value to a list.
        prob_distribution.append(prob)

    for i in range(n - k + 1):
        prob_distribution[i] = prob_distribution[i]/prob_sum
        #print(i, prob_distribution[i])

        """Now we will take a break in the middle of this function to move to part 2 of the function.
        Here we will call function 'roll_weighted_dice()' that will accept this distribution list of
        probabilities. """

    # Call the function that will return the kmer index required to grab the right kmer and plug
    # it back into our kmer_list.
    selection_number = roll_weighted_dice(prob_distribution)

    
    selected_kmer = kmer_list[selection_number]

    return selected_kmer
        
        

def roll_weighted_dice(probabilities):

    """ The function will add a bit of randomness into into picking a kmer from this list,
    but will still bias the randomness towards kmers with higher probabilities.  Without this added
    randomness, the function would simply output the kmer with the highest prob score."""

    # Create a list of sides based on the number of probabilities
    sides = numpy.arange(0, len(probabilities))
        
    # Use numpy.random.choice to simulate the die roll
    result = numpy.random.choice(sides, p = probabilities)

    return result



def score(motifs, k, t):
    """ Calculate the number of mismatches between a list of kmers and the list's consensus kmer"""
    
    consensus = ""

    # Rather than importing a profile from the function, temporarily construct a profile matrix
    # from the motifs list and use that to generate the consensus string.
    profile = form_profile(motifs, k)

    # iterate over each position in the matrix profile. 
    for i in range(k):
        # iterate over the ACGT keys in the ith row.
        for key in profile[i]:  
            # Grab the max probability corresponding to a given base key.
            m = max(profile[i].values())
            # If a profile key in the i-th row is == to max, add it to the consensus strand.
            if profile[i][key] == m:
                consensus += key
                break # Break so that we don't add >1 base in cases of equal probabilities for diff. keys. 
    #print('consensus', consensus)

    
    # Baseline score = 0. 
    score = 0
    # iterate over kmers in list. 
    for motif in motifs:
        #print('motif', motif)
        for i in range(k):
            
            # Count up for ever mismatch agains the consensus strand.
            if motif[i] != consensus[i]:
                score += 1
            #print(motif[i], consensus[i], score)
        #print(motif, score)
    return score





