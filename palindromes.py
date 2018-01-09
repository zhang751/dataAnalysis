def is_palindrome(s1):
    ''' (str) -> bool
    
    Return True iff s1 is a palindrome, which means that s reads the same 
    forwards and backwards.
    
    Precondition: s1.islower() == True
    
    >>> is_palindrome('abcdefg')
    False
    >>> is_palindrome('abcdcba')
    True
    >>> is_palindrome('')
    True
    '''
    
    reverse_s1 = ''
    for char in s1:
        reverse_s1 = char + reverse_s1
    return s1 == reverse_s1


def is_palindromic_phrase(s2):
    ''' (str) -> bool
    
    Return True iff s2 is a palindrome ignoring all non-alphabetic characters, 
    and that uppercases are considered to be equivalent to lowercases.
    
    >>> is_palindromic_phrase('acw23rgsHJG')
    False
    >>> is_palindromic_phrase('3aCk23d4d1Kc4a')
    True
    >>> is_palindromic_phrase('')
    True
    '''
    
    x = ''
    # Create a new strand that does not contain any capital letters or digits.
    for i in range(len(s2)):
        if s2[i].islower():
            x = x + s2[i]
        elif s2[i].isupper():
            x = x + s2[i].lower()
        else:
            x = x
    return is_palindrome(x)
    
    
def get_odd_palindrome_at(s3, center_index):
    ''' (str, int) -> str
    
    Precondition: s3.islower() == True and s3[center_index] is valid.
    
    Return the longest odd-length palindrome in s3 centered at center_index.
    
    >>> get_odd_palindrome_at('aabbcbbaa', 4)
    'aabccbbaa'
    >>> get_odd_palindrome_at('abacdfex', 0)
    'a'
    >>> get_odd_palindrome_at('xyxcmofem', 1)
    'xyx'
    '''
    
    i = center_index
    n = 0
    # Return s3[i] if i is the index of the first character.
    if i == 0:
        return s3[0]
    else:
    # Check each character at both left and right side starting from 
    # the center_index.
        for k in range(min(len(s3) - i, i + 1)):
            if s3[i - k] == s3[i + k]:
                n = n + 1
            else:
                k = min(len(s3) - i, i)
        return s3[i - n + 1:i + n]