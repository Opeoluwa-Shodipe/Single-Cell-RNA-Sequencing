"""
Hackbio Internship
Team: Glycine

# Task: Calculate Hamming Distance
# Author: Opeoluwa Shodipe
# Github: https://github.com/Opeoluwa-Shodipe/Single-Cell-RNA-Sequencing/tree/main/Stage%20one
# Linkedin: https://www.linkedin.com/in/sopeoluwa/
"""

def hamming_distance(str1, str2, ignore_case=True):
    """
    Compute the Hamming distance between two strings of equal or unequal length.

    Args:
        str1 (str): First string.
        str2 (str): Second string.
        ignore_case (bool): If True, compare case-insensitively.

    Returns:
        int: Count of differing positions.
    """
    if ignore_case:
        str1, str2 = str1.lower(), str2.lower()
    max_len = max(len(str1), len(str2))
    padded1 = str1.ljust(max_len)
    padded2 = str2.ljust(max_len)
    return sum(ch1 != ch2 for ch1, ch2 in zip(padded1, padded2))

distance = hamming_distance("Eadencre8ives", "eadencre8ives")
print(distance)  # Result: 0 (if ignore_case=True)

