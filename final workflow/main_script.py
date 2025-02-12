import pandas as pd
from multiprocessing import Pool, cpu_count
from alignment_utils import generate_sorted_substrings, align_sequence
from using_meme import meme_analysis
from fasta_converter import fasta_format_converter
from seed_finder import get_top_sequences, seed_meme_analysis, read_html_pwm, calculate_consensus

def read_sequences_from_fasta_file(file_path, col_index):
    with open(file_path, 'r') as file:
        next(file)
        sequences = []
        scores = []
        for line in file:
            if line.startswith(">"):
                continue
            else:   
                if col_index == 0:
                    sequences.append(line.strip())
                else:
                    sequences.append(line.strip())
    return sequences 

def read_sequences_from_file(file_path, col_index):
    with open(file_path, 'r') as file:
        next(file)
        sequences = []
        scores = []
        for line in file:
            columns = line.split()
            if col_index == 0:
                sequences.append(columns[col_index])
            else:
                sequences.append(columns[col_index][::-1])
            scores.append(float(columns[2]))
    return sequences, scores

def align_sequences_parallel(sequences, reference_sequence):
    substrings = generate_sorted_substrings(reference_sequence)
    ref_index = 0
    tasks = [(seq, reference_sequence, substrings, ref_index) for seq in sequences]
    
    with Pool(processes=cpu_count()) as pool:
        dfs = pool.map(align_sequence, tasks)

    # Combine all partial DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(sequences, 1):
            file.write(f">sequence_{i}\n{sequence}\n")
       

if __name__ == "__main__":
    for TF in ['Jun_Fos']:
        for col_index in [0,1]:
            input_file_path = f'MonomerBinding/{TF}/{TF}.txt'
            sequences, scores = read_sequences_from_file(input_file_path, col_index)
            
            top_sequences = get_top_sequences(sequences, scores, top_n=20)        
            
            output_file_path = f'MonomerBinding/{TF}/col_{col_index}/{TF}_consensus.txt'
            output_fasta_path = f'MonomerBinding/{TF}/col_{col_index}/{TF}_consensus.fasta'

            output_fasta_file = f'MonomerBinding/{TF}/col_{col_index}/{TF}_top_20_sequences.fasta'
            write_sequences_to_fasta(top_sequences, output_fasta_file)
            
            print("****************************************************************")
            seed_meme_analysis(output_fasta_file, col_index, 20)
            print("****************************************************************")
            html_file = f'MonomerBinding/{TF}/col_{col_index}/Meme_of_top_20_Seeds/meme.html'
            pwm_section = read_html_pwm(html_file)
            reference_sequence = calculate_consensus(pwm_section)
            
            length_of_sequence = len(sequences[0])

            print(reference_sequence, ' ', col_index)
            if not reference_sequence:
                raise ValueError(f"No sequence containing '{reference_sequence}' was found in the file.")

            df_aligned = align_sequences_parallel(sequences, reference_sequence)

            columns_to_extract = [col for col in range(length_of_sequence) if col in df_aligned.columns]
            extracted_data = df_aligned[columns_to_extract]
            extracted_data = extracted_data.fillna('-')

            formatted_data = extracted_data.apply(lambda row: ''.join(row.astype(str)), axis=1)

            print("**** Output File ***")
            with open(output_file_path, 'w') as file:
                file.write(reference_sequence + '\n')
                for line in formatted_data:
                    file.write(line + '\n')
                    
            # sequence_logo_generator(TF,output_file_path,col_index)
            fasta_format_converter(output_file_path, output_fasta_path)
            meme_analysis(output_fasta_path, length_of_sequence)