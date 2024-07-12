import sys

# NOTE: THIS CODE SUCCESSFULLY PASSES ALL TEST FOR WEEK 1 CHAPTER 8.4.
    
# Insert your neighbors function here, along with any subroutines you need
def neighbors(s, d):

    if d == 0:
        return [s]    # Base case
    
    bases = ['A', 'C', 'T', 'G']    # out list of bases
    neighborhood = set([s])        # Empty list except for original sequence

    for i in range(len(s)):    #iterate through s
        for base in bases:
            if s[i] != base:
                
                neighbor = s[:i] + base + s[i+1:]
                
                other_neighbors = neighbors(neighbor, d - 1)
                
                neighborhood.update(other_neighbors)

    return list(neighborhood)




def frequent_words_with_mismatches(text, k, d):
    """Find the most frequent k-mers with up to d mismatches in a text."""

    n = len(text)
    patterns = set()   # empty list to store patterns? Do I need this?
    freqMap = {}  #dict for patterns corresponding to frequency in text.

    mismatches = []

    

    # Call neighbors
    for i in range(n-k+1):
        
        k_mer = text[i:i+k]     # grab k chars at a time. 
        #print('actual k-mer', pattern)  # Test print

        for j in neighbors(k_mer, d):

            if j in patterns:
                freqMap[j] += 1
            else:
                freqMap[j] = 1

            patterns.update([j])    # Update the pattern into the set of neighbors
            # Don't forget to add brackets [] around j otherwise the set will be individ. chars.
            # Such as A, T, C, G only, not whole patterns.

    m = max(freqMap.values())   # Get the mismatch maximum
    
    for key in freqMap:
        if freqMap[key] == m:

            mismatches.append(key)   # Create a new list of the final patterns. 
            

    return mismatches


text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
d = 1
k = 4
print(frequent_words_with_mismatches(text, k, d))

