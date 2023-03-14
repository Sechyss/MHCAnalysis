import pandas as pd


def transformmhclist(filename):
    df = pd.read_table(filename, delim_whitespace=True, header=None)
    listmhc = list(df.iloc[:, 1])
    listlength = list(df.iloc[:, 2])

    return listmhc, listlength


def generate_sequences(reference_seq, snps):
    """Generate all potential sequences based on a reference sequence and SNPs"""
    sequences = [reference_seq]
    for snp in snps:
        new_sequences = []
        for seq in sequences:
            # Apply SNP to the sequence
            new_seq = list(seq)
            new_seq[snp['pos']] = snp['alt']
            new_seq = ''.join(new_seq)
            new_sequences.append(new_seq)
        sequences.extend(new_sequences)
    return sequences


def generate_combinations(ref_seq, snp_dict):
    """
    Generate all possible combinations of the reference sequence and SNPs dictionary.

    Parameters:
    ref_seq (pd.Series): Pandas Series representing the reference DNA sequence.
    snp_dict (dict): Dictionary with SNPs using the position as keys and a list of nucleotides as the value.

    Returns:
    list: List of all possible combinations as strings.
    """
    # Initialize the list of all possible combinations with the reference sequence
    combinations = ["".join(ref_seq)]

    # Loop through the SNP dictionary and generate all possible combinations
    for pos, nts in snp_dict.items():
        # Check that the SNP position is within the range of valid indices for the reference sequence
        if pos < 0 or pos >= len(ref_seq):
            raise ValueError(f"Invalid SNP position: {pos}")
        # Make a copy of the list of combinations
        new_combinations = combinations.copy()
        # Loop through the list of combinations and generate a new combination for each SNP
        for nt in nts:
            for y in range(len(combinations)):
                # Replace the nucleotide at the SNP position with the SNP nucleotide
                new_seq = list(combinations[y])
                new_seq[pos] = nt
                # Add the new combination to the list of new combinations
                new_combinations.append("".join(new_seq))
        # Update the list of combinations with the new combinations
        combinations = new_combinations

    return combinations
