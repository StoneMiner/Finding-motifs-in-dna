import sys
import random
import subprocess

""" this program attempts to find DosR transcription factor binding sites 
in tuberculosis-causing bacteria. """

# Insert your randomized_motif_search function here, along with any subroutines you need
def read_input_file(input_file_path, t):
    """ Funciton reads a variety of input txt files to pass into this script. """
    with open(input_file_path, 'r') as file:
        # Grab the dna text" 
        content = file.read()
        content = content.replace('\n', '')

    #make an empty list for copies of this dna
    dna = []
    
    for i in range(5):
        
        dna.append(content)
    
        
    return dna


def find_dosr_motif(dna):
    """ For now, call a gibbsSampler program to generate some possible motifs. """

    #k = random.randint(10, 30)
    k = 13
    
    n = 100
    
    from modGibbsSampler import gibbs_sampler

    best_kmers, best_score = gibbs_sampler(dna, k, t, n)

    print(best_kmers, best_score)

t = 5
dna = read_input_file('input.txt', t)

find_dosr_motif(dna)
