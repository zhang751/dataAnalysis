def is_base_pair(base1, base2):
    ''' (str, str) -> bool
    
    Precondition: base1 and base2 can only be 'A' or 'T' or 'G' or 'C'
    
    Return True iff base1 and base2 can form a DNA base pair.
    
    >>> is_base_pair('A', 'T')
    True
    >>> is_base_pair('A', 'C')
    False
    '''
    
    x = base1 + base2
    if 'A' in x and 'T' in x:
        return True
    elif 'C' in x and 'G' in x:
        return True
    else:
        return False
    
    
def is_dna(strand1, strand2):
    ''' (str, str) -> bool
    
    Precondition: len(strand1) == len(strand2) and 
    strand1 and strand2 can only contain 'A', 'T', 'C', 'G'
    
    Return True iff strand1 and strand2 can properly form a DNA molecule.
    
    >>> is_dna('ACTG', 'TGAC')
    True
    >>> is_dna('ACTTAG', 'TCCCGC')
    False
    '''
    
    # Test if the corresponding indices in strand1 and strand2 form base pairs.
    n = 0
    for i in range(len(strand1)):
        if is_base_pair(strand1[i], strand2[i]):
            n = n + 1
    return n == len(strand1)


def is_dna_palindrome(strand3, strand4):
    ''' (str, str) -> bool
    
    Precondition: is_dna(strand3, strand4) == True
    
    Return True iff strand3 and strand4 can form a DNA palindrome
    
    >>> is_dna_palindrome('ATCTA', 'TAGAT')
    False
    >>> is_dna_palindrome('TGCAGT', 'ACGTCA')
    False
    >>> is_dna_palindrome('CCC', 'CCC')
    False
    '''
    
    # Reverse strand3 to see if it is the the same as strand4.
    reverse_strand3 = ''
    if is_dna(strand3, strand4):
        for char in strand3:
            reverse_strand3 = char + reverse_strand3
        return reverse_strand3 == strand4
    else:
        return False
    
    
def restriction_sites(dna_strand1, recog_seq):
    '''  (str, str) -> list of int
    
    Return a list of all indices in dna_strand that where recog_seq appears.
    
    >>> restriction_sites('GCCTAGCCAGCTGATTCGATTCAGCTGA', 'CCGCTA')
    []
    >>> restriction_sites('GGGGGGG', 'GG')
    [0, 2, 4]
    >>> restriction_sites('GCATCATAAA', 'GC')
    [0]
    '''
    
    sites = []
    x = 0    
    if recog_seq in dna_strand1:
        while x < (len(dna_strand1) - len(recog_seq)):
    # To avoid the dead loop if .find method evaluates to -1.
            if dna_strand1.find(recog_seq, x) == -1:
                return sites
            else:
                sites.append(dna_strand1.find(recog_seq, x))
                x = dna_strand1.find(recog_seq, x) + len(recog_seq)
        return sites
    # If there's no recog_seq in dna_strand1, return the empty list.
    else:
        return []


def match_enzymes(dna_strand2, res_ems1, rec_seqs1):
    ''' (str, list of str, list of str) -> list of two-item 
    [str, list of int] list
    
    Return a two-item list where the first item representing emzymes
    in res_ems1, and the second item representing the indices of
    dna_strand2 where the corresponding rec_seqs1 are found.
    
    >>> match_enzymes('GCTAGCTGCGGATCCGTTA', ['BamHI'], ['GGATCC'])
    [['BamHI', [9]]]
    >>> match_enzymes('AGCTTCGATCGAAGTAGCT', ['TaqI', 'AluI'], ['TCGA', 'AGCT'])
    [['TaqI', [4, 8]], ['AluI', [0, 15]]]
    '''
    
    emzs = []
    i = 0
    while i < len(res_ems1):
    # Call restriction_sites function to find the indices, and then
    # append the list of indices of each emzyme to the returning list.
        emzs.append([res_ems1[i], restriction_sites(dna_strand2, rec_seqs1[i])])
        i = i + 1
    return emzs


def one_cutters(dna_strand3, res_ems2, rec_seqs2):
    ''' (str, list of str, list of str) -> list of two-item [str, int] list
    
    Return a two-item list where the first item is the name of the enzyme in 
    res_ems2 and the second item is the index of dna_strand3 where the 
    corresponding rec_seqs1 are found. Return an empty string
    
    >>> one_cutters('GCTAGCTGCGGATCCGTTA', ['BamHI'], ['GGATCC'])
    [['BamHI', 9]]
    >>> one_cutters('AGCTTCGATCGAAGTAGCT', ['TaqI', 'AluI'], ['TCGA', 'AGCT'])
    []
    '''
    
    all_onecutters = []
    i = 0
    while i < len(rec_seqs2):
    # Check if the rec_seq only appears once in dna_strand3.
        if dna_strand3.count(rec_seqs2[i]) > 1:
            i = i + 1
        elif dna_strand3.count(rec_seqs2[i]) == 0: 
            i = i + 1
        else:
    # Collect all indices of one cutters in dna_strand3 into a list
            all_onecutters.append([res_ems2[i], dna_strand3.find(rec_seqs2[i])])
            i = i + 1
    return all_onecutters


def correct_mutations(mutated_dna, clean, res_ems3, rec_seqs3):
    ''' (list of str, str, list of str, list of str) -> NoneType
    
    Precondition: the clean strand only contains 1 one-cutter from the res_ems3
    and rec_seqs3
    
    The function modifies mutated_dna that shares exactly 1 one-cutter with
    the clean strand by replacing all bases starting at the one-cutter in the 
    mutated_dna strand with all bases starting at the one-cutter in the clean 
    strand, up to and including the end of the strand.
    
    >>> mutated_dna = ['AACTTGCTAGCTAAC', 'TTCCAACCTT']
    >>> clean = 'CCGAGCTAAGC'
    >>> res_ems3 = ['HgaI', 'AluI']
    >>> rec_seqs3 = ['GACGC', 'AGCT']
    >>> correct_mutations(mutated_dna, clean, res_ems3, rec_seqs3)
    >>> mutated_dna
    ['AACTTGCTAGCTAAGC', 'TTCCAACCTT']
    
    >>> mutated_dna = ['AGATCAACGGA', 'CCATTGCAT']
    >>> clean = 'CCAGATCAACT'
    >>> res_ems3 = ['Sau3A', 'EcoRI']
    >>> rec_seqs3 = ['GATC', 'GAATTC']
    >>> correct_mutations(mutated_dna, clean, res_ems3, rec_seqs3)
    >>> mutated_dna
    ['AGATCAACT', 'CCATTGCAT']
    '''
    
    # Modify each string in mutated_dna one by one.
    for i in range(len(mutated_dna)):
        x = ''
        strand = mutated_dna[i]
        n = 0
        while n < len(rec_seqs3):
    # Determine if the sequence appears in the clean strand and that the 
    # sequence only appears once in mutated_dna
            if mutated_dna[i].count(rec_seqs3[n]) == 1 and \
               rec_seqs3[n] in clean:
                a = strand.find(rec_seqs3[n])
                b = clean.find(rec_seqs3[n])
                x = strand[:a] + clean[b:]
                mutated_dna[i] = x
            else:
                mutated_dna[i] = mutated_dna[i]
            n = n + 1   
            
if __name__ == '__main__':
    import doctest
    doctest.testmod()
