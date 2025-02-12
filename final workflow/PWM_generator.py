import numpy as np

def count_bases(filename):
    with open(filename, 'r') as file:
        sequences = [line.strip().upper() for line in file.readlines()]
    
    num_sequences = len(sequences)
    sequence_length = 10

    count_matrix = np.zeros((5, sequence_length), dtype=int)

    base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}

    for sequence in sequences:
        for i, base in enumerate(sequence):
            if base in base_index:
                count_matrix[base_index[base]][i] += 1

    probability_matrix = count_matrix / num_sequences

    return probability_matrix

def save_matrix_to_file(matrix, output_filename):
    base_labels = ['A', 'C', 'G', 'T', '-']  # Base labels

    with open(output_filename, 'w') as f:
        for i, row in enumerate(matrix):
            f.write(base_labels[i] + '\t' + '\t'.join(f'{value:.2f}' for value in row) + '\n')

for TF in ['FOS_CEBPE_20N']:
    for col_index in [0,1]:
        input_filename = f'/Users/bristi/Desktop/Design Project/Working-with-TF/DimerBinding/{TF}_col_{col_index}.txt'
        output_filename = f'/Users/bristi/Desktop/Design Project/Working-with-TF/DimerBinding/PWM_{TF}_col_{col_index}.txt'

        prob_matrix = count_bases(input_filename)
        save_matrix_to_file(prob_matrix, output_filename)

        print(f"Probability matrix saved to {output_filename}")