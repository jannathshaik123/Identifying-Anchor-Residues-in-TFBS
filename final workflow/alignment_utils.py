import pandas as pd
from tqdm import tqdm

def generate_sorted_substrings(s):
    substrings = [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    return unique_substrings

def align_sequence(sequence_info):
    sequence, reference_sequence, substrings, ref_index = sequence_info
    max_length = len(sequence[0])
    lower_bound = 1 - max_length
    upper_bound = 2 * max_length - 1
    columns = [i for i in range(lower_bound, upper_bound)]
    df = pd.DataFrame(columns=columns)

    aligned = False
    if reference_sequence in sequence:
        idx = sequence.find(reference_sequence)
        aligned = True
        entire = True
    else:
        prev_best_number_of_matches = 0
        best_substring_index = -1

        # Add tqdm for progress tracking in the inner loop
        for i in tqdm(range(len(substrings)), desc="Finding best substring match"):
            total_number_of_matches = 0
            
            if substrings[i] in sequence:
                current_substring = substrings[i]
                idx = sequence.find(current_substring)
                substring_index_in_seed = reference_sequence.find(current_substring)
                offset = substring_index_in_seed - idx

                if offset < 0:
                    aligned_seed = reference_sequence[:len(reference_sequence) - abs(offset)]
                    aligned_sequence = sequence[abs(offset):]
                else:
                    aligned_seed = reference_sequence[abs(offset):]
                    aligned_sequence = sequence[:len(sequence) - abs(offset)]

                min_length = min(len(aligned_seed), len(aligned_sequence))
                total_number_of_matches = sum(1 for i in range(min_length) if aligned_seed[i] == aligned_sequence[i])

                if prev_best_number_of_matches < total_number_of_matches:
                    best_substring_index = i
                    prev_best_number_of_matches = total_number_of_matches
        
        best_substring = substrings[best_substring_index]
        idx = sequence.find(best_substring)
        substring_index_in_seed = reference_sequence.find(best_substring)
        aligned = True
        entire = False

    if aligned:
        if entire:
            offset = ref_index - idx
        else:
            offset = substring_index_in_seed - idx
        indices = range(offset, offset + len(sequence))
        df = pd.DataFrame([list(sequence)], columns=indices)

    return df
