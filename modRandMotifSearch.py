import sys
import random

sys.setrecursionlimit(2000)
# Please do not remove package declarations because these are used by the autograder.

# Insert your randomized_motif_search function here, along with any subroutines you need
def read_input_file(input_file_path):
    """ Funciton reads a variety of input txt files to pass into this script. """
    with open(input_file_path, 'r') as file:
        # Read the first line for k and t
        first_line = file.readline().strip()

        k, t = map(int, first_line.split())

        # Grab the dna text and create individual list elements based on spaces. 
        dna = []
        for line in file.readlines():
            line_dna = line.strip().split()
            dna.extend(line_dna)  # Flatten the list of lists into a single list

    return dna, k, t

def randomized_motif_search(dna: list[str], k: int, t: int) -> list[str]:
    """Implements the RandomizedMotifSearch algorithm with pseudocounts.

    Step 1. Use random number generator to pick random kmers from each
    string of dna.
    2. collect these kmers into a list.
    3. Create a profile probability matrix from this list.
    4. Use the profile matrix to assess the most probable kmer from each
        dna string, forming a new list of kmers.
    5. Score your random kmers and score your matrix-selected kmers.
    6. Restart the algorithm, going all the way back to generating another random set of
        motifs, and using these random motifs to build a profile matrix that pulls out a
        set of motifs.  You will score them again and continue to see if we can get
        ever lower scores.
    7. Do this 1000 times.

    Update. The thing that this code is currently doing: I start of with an initial set of
    random kmers as best motifs.
    I make a matrix profile from this list, then use it to find_most_probable_kmer.
    I score the matrix selected kmers and the random ones, and see which is better.
    I better one gets saved in best_motifs.

    What next? is more randomness supposed to be introduced into the equation? Or is the random
    first set of kmers enough to start with, and then after that we can dispense with randomness
    and only use iteratively better and better kmers and kmer profile matricies?

    Different motifs will be scored differently by different profile matricies.  This is problematic
    if the some matricies are built on 'bad' motif lists.  Maybe when we have completed our initial
    score to get best_kmers, we can build our new profile matrix on best_kmers, then restart our
    randomness loop.

    Here's what we will do with the new random kmers. We won't automatically make the new random kmers
    into a profile.  We will score them based on our existing profile matrix (at this point this is profile
    matrix #2).  If the random motifs score better than our profile, we can make a new profile.
    If not, we can keep generating random kmers and scoring them until a better one comes.


    Here's another thing.  Maybe I've been thinking about this problem all wrong.  Maybe continuously improving
    the score of 'most_probable_kmers' isn't ideal, because the scoring criteria will always be provided
    only by the first matrix profile, or a derivitive of it, which is itself random.
    Perhaps what we really want is to see which set of random kmers are capable of forming their own matrix profile
    and their own most_probable_mers that have the lowest overall score. In other words.  We start the loop before
    making randoms, and we don't restart the loop until after most_probable_kmers have been selected.
    
    """
    loop_counter = 1000

    lowest_score = 1000

    best_motifs = []

    while loop_counter > 0:
    #while lowest_score > 5:
        loop_counter -= 1
        # Call randomizer and get a list of kmers
        randoms = get_random_kmers(dna, k)

        # Get our profile matrix from the random motifs
        matrix_profile = form_profile(randoms, k)
        
        # Updated list calls the improve_motifs function, passing in the dna list, and the current kmer list.
        selected_motifs = []
        for i in range(t):
            dna_string = dna[i]
            selected_motifs.append(find_most_probable_kmer(dna_string, matrix_profile, k))


        # At this point in the function we have 1) A profile matrix formed from random kmers,
        # and 2) a collection of matrix-selected kmers, one from each dna string (= t kmers).

        # Now we can proceed to score the random kmers and the selected kmers.
        # Score both the lists
        score_randoms = score(randoms, k, t)
        score_selected = score(selected_motifs, k, t)

        # Best motifs will always be the one with lower score.
        if score_selected <= score_randoms:
            
            if score_selected < lowest_score:
                lowest_score = score_selected
                best_motifs = selected_motifs

        #print(score_selected, selected_motifs)

    # Call the comparison to output function
    #print(compare_output(best_motifs, dna))

    return lowest_score, best_motifs
        


def get_random_kmers(dna: list[str], k: int) -> list[str]:
    """ This function randomly selects kmers in each string of dna"""
    
    # length of a given dna string
    n = len(dna[0])
    kmer_count = n - k - 1

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

            
def find_most_probable_kmer(dna, profile, k):
    """ Based on a matrix of base probabilities and an input dna string, this function will
    go over every kmer in the dna string and find the one with the closest matching probability
    to the matrix"""
    # Initialize max_prob with bencmark low value, a negative.
    max_prob = -1
    # Default/first most probable kmer will simply be first kmer in the dna string, to start off. 
    most_probable_kmer = dna[:k]

    n = len(dna)
    # iterate over the dna string and pull out kmers. 
    for i in range(n - k + 1):
        kmer = dna[i:i + k] # each kmer.
        
        prob = 1    # Default probability will be 1 so we can multiply.
        # for each position in kmer, and each character in kmer (each base),  multiply the probability
        # against the corresponding probability in our matrix profile for that position.
        # store the new probability value in 'prob', our new benchmark. 
        for j, char in enumerate(kmer):
            
            prob *= profile[j][char]
        # Iteratively sort our highest probability and store in max_prob.
        if prob > max_prob:
            max_prob = prob
            # The kmer that produces the highest probability is the most_probable. 
            most_probable_kmer = kmer
    
    #print('most_probable_kmer', most_probable_kmer)        
    return most_probable_kmer
        


def score(motifs, k, t):
    """ Calculate the number of mismatches between a list of kmers and the list's consensus kmer"""
    consensus = ""

    # Form the motifs list and use that to generate the consensus string.
    profile = form_profile(motifs, k)

    # iterate over each position in the matrix profile. 
    for i in range(k):
        # iterate over the ACGT keys in the ith row.
        #print(i, profile[i])                # Test print
        for key in profile[i]:  
            # Grab the max probability corresponding to a given base key.
            m = max(profile[i].values())
            # If a profile key in the i-th row is == to max, add it to the consensus strand.
            if profile[i][key] == m:
                consensus += key
                break # Break so that we don't add >1 base in cases of equal probabilities for diff. keys. 

    
    # Baseline score = 0. 
    score = 0
    # iterate over kmers in list. 
    for motif in motifs:
        #print('motif', motif)
        for i in range(k):
            # Count up for ever mismatch agains the consensus strand.
            if motif[i] != consensus[i]:
                score += 1
    return score


def read_output_file(output_file_path):
    """ this function reads output files and grades the assignment """
    with open(output_file_path, 'r') as file:
        #grab the kmer lists on the first line.
        answer_list = []
        for line in file.readlines():
            kmer = line.strip().split()
            answer_list.extend(kmer)

    return answer_list

def compare_output(output: list[str], dna) -> str:
    """Grades the output of the code with the problem set solution"""

        # Grab the output file list kmers for the correct answer, to compare. 
    solution = read_output_file(output_file_path)
    
    #Base case
    if output == None:
        return
    
    # Script to compare our functions output kmers with the answer from the text file.
    match_counter = 0
    for i in range(len(output)):
        print(output[i], solution[i])
        if output[i] == solution[i]:
            match_counter += 1
            print(match_counter)

            
    if match_counter == len(dna):
        return 'correct output'
    else:
        return 'wrong output'

            
input_file_path = 'RandomizedMotifSearch/inputs/input_1.txt'  #File path for debug data. Replace asneeded.
output_file_path = 'RandomizedMotifSearch/outputs/output_1.txt'
dna, k, t = read_input_file(input_file_path)  # Extract dna, k, and t.
print(randomized_motif_search(dna, k, t))
