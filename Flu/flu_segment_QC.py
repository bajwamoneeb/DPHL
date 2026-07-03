import os
import subprocess
import datetime
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# Expected segments for A-type and B-type influenza
A_TYPE_SEGMENTS = ['A_HA_H1', 'A_HA_H2', 'A_HA_H3', 'A_MP', 'A_NA_N1', 'A_NA_N2', 'A_NA_N3', 'A_NP', 'A_NS', 'A_PA', 'A_PB1', 'A_PB2']
B_TYPE_SEGMENTS = ['B_HA', 'B_MP', 'B_NA', 'B_NP', 'B_NS', 'B_PA', 'B_PB1', 'B_PB2']

# Desired output column order
DESIRED_COLUMNS = ['Seq_Id (HA)', 'Seq_Id (NA)', 'Seq_Id (PB1)', 'Seq_Id (PB2)',
                   'Seq_Id (PA)', 'Seq_Id (MP)', 'Seq_Id (NS)', 'Seq_Id (NP)',
                   'Seq_Id (HE)', 'Seq_Id (P3)']

# Function to adjust record IDs by removing prefixes before flu type
def adjust_record_id(record_id):
    # If the record ID contains an underscore and the flu type identifier (A or B)
    if '_' in record_id:
        parts = record_id.split('_')
        # Find the index where 'A' or 'B' appears
        for i, part in enumerate(parts):
            if part in ['A', 'B']:
                # Return the rest of the ID starting from 'A' or 'B'
                return '_'.join(parts[i:])
    return record_id  # Return the original if no adjustment is needed

# Function to run BLAST for a single segment and check PID
def run_blast_for_segment(record, prefix, low_pid_file, blast_results_file, db="nt"):
    try:
        result = subprocess.run([
            'blastn', '-query', '-', '-db', db, '-remote', '-outfmt', '6 qseqid sseqid pident',
            '-max_target_seqs', '1', '-evalue', '1e-5'
        ], input=str(record.seq), text=True, capture_output=True)

        # Log the BLAST results to the separate file
        with open(blast_results_file, 'a') as brf:
            brf.write(f'{datetime.datetime.now()} - BLAST Results for {record.id} in {prefix}:\n')
            brf.write(result.stdout + '\n')

        pid_values = [float(line.split()[2]) for line in result.stdout.strip().split('\n') if line]
        if pid_values and max(pid_values) < 99.0:
            with open(low_pid_file, 'a') as lpf:
                lpf.write(f'{datetime.datetime.now()} - Segment {record.id} for {prefix} has a top match below 99%: {max(pid_values)}%\n')
        return record.id, result.stdout.strip()
    except Exception as e:
        return f'Error running BLAST for {record.id}: {str(e)}\n', None

# Function to process each fasta file
def process_fasta_file(fasta_file, low_pid_file, blast_results_file, log_file):
    prefix = os.path.basename(fasta_file).split('.')[0].split('-')[0]
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    errors = []
    results = []
    with ProcessPoolExecutor() as executor:
        # Pass all required arguments to run_blast_for_segment
        futures = [
            executor.submit(run_blast_for_segment, record, prefix, low_pid_file, blast_results_file)
            for record in records
        ]
        for future in as_completed(futures):
            error, result = future.result()
            if error:
                errors.append(error)
            if result:
                results.append(result)
    return errors, results

# Function to convert fasta file to dataframe
def convert_fasta_to_df(fasta_file):
    # Adjust headers when reading the records
    fasta_sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        adjusted_id = adjust_record_id(record.id)
        record.id = adjusted_id
        record.description = ''
        fasta_sequences.append(record)
    headers = [record.id for record in fasta_sequences]
    sequences = [str(record.seq) for record in fasta_sequences]
    if len(headers) != len(sequences):
        raise ValueError("Headers and sequences lists are of different lengths")
    file_prefix = os.path.basename(fasta_file).split('.')[0].split('-')[0]
    df = pd.DataFrame({'headers': headers, 'sequences': sequences})
    transposed_df = df.set_index('headers').T
    transposed_df['File_Prefix'] = file_prefix
    return transposed_df

def process_low_pid_log(log_file):
    """Parse the low PID log to identify segments."""
    low_pid_segments = set()
    try:
        with open(log_file, 'r') as lf:
            for line in lf:
                if 'below 99%' in line:
                    parts = line.strip().split()
                    prefix = parts[6]
                    segment = parts[4]
                    segment_id = f"{prefix}-{segment}"
                    low_pid_segments.add(segment_id)
        return low_pid_segments
    except FileNotFoundError:
        return set()

# Function to compile sequences, create key file, and filter based on low PID
def compile_and_filter_sequences(input_dir, low_pid_file, output_fasta_path, output_tsv_path, key_file_path, counter_file):
    df_list = []
    seq_records = []
    id_key = []

    if os.path.exists(counter_file):
        with open(counter_file, 'r') as cf:
            record_counter = int(cf.read().strip())
    else:
        record_counter = 1

    prefix_counter = {}

    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith('.fasta'):
            try:
                df = convert_fasta_to_df(os.path.join(input_dir, fasta_file))
                df_list.append(df)
                # Adjust headers when reading the records
                records = []
                for record in SeqIO.parse(os.path.join(input_dir, fasta_file), 'fasta'):
                    adjusted_id = adjust_record_id(record.id)
                    record.id = adjusted_id
                    record.description = ''
                    records.append(record)
                prefix = os.path.basename(fasta_file).split('.')[0].split('-')[0]
                for record in records:
                    if prefix not in prefix_counter:
                        prefix_counter[prefix] = record_counter
                        record_counter += 1
                    new_id = f"DE-DHSS-{prefix_counter[prefix]}_{record.id}"
                    seq_records.append(SeqRecord(record.seq, id=new_id, description=""))
                    id_key.append((f"{prefix}-{record.id}", new_id))
            except Exception as e:
                print(f"Error processing {fasta_file}: {e}")

    if df_list:
        combined_df = pd.concat(df_list, ignore_index=True)

    if seq_records:
        with open(output_fasta_path, 'w') as fasta_output:
            SeqIO.write(seq_records, fasta_output, "fasta")

        with open(key_file_path, 'w') as key_output:
            pd.DataFrame(id_key, columns=['Original_ID', 'New_ID']).to_csv(key_output, sep='\t', index=False)

    # Ensure `process_low_pid_log` handles empty or missing files gracefully
    low_pid_segments = process_low_pid_log(low_pid_file)

    id_key_df = pd.read_csv(key_file_path, sep='\t')
    id_key_dict = dict(zip(id_key_df['Original_ID'], id_key_df['New_ID']))

    filtered_records = []
    for record in seq_records:
        original_id = next((k for k, v in id_key_dict.items() if v == record.id), None)
        if original_id and original_id not in low_pid_segments:
            filtered_records.append(record)

    if filtered_records:
        with open(output_fasta_path, 'w') as fasta_output:
            for record in filtered_records:
                fasta_output.write(f">{record.id}\n{str(record.seq)}\n")
    else:
        print("No valid data to write to FASTA.")

    if filtered_records:
        filtered_headers = [record.id for record in filtered_records]
        filtered_id_key_df = id_key_df[id_key_df['New_ID'].isin(filtered_headers)].copy()
        if not filtered_id_key_df.empty:
            filtered_id_key_df['File_Prefix'] = filtered_id_key_df['Original_ID'].str.split('-').str[0].str.split('.').str[0]
            filtered_id_key_df['Segment'] = filtered_id_key_df['Original_ID'].str.split('-').str[1]
            df_pivot = filtered_id_key_df.pivot(index='File_Prefix', columns='Segment', values='New_ID').reset_index()

            # Ensure the filtered_id_key_df has the correct column headers and order for both A-type and B-type
            df_pivot = df_pivot.rename(columns=lambda x: 'Seq_Id (' + x.split('_')[1] + ')' if x != 'File_Prefix' else x)

            # Reorder the columns according to DESIRED_COLUMNS
            column_order = ['File_Prefix'] + DESIRED_COLUMNS
            for col in column_order:
                if col not in df_pivot.columns:
                    df_pivot[col] = None
            df_pivot = df_pivot[column_order]

            df_pivot.to_csv(output_tsv_path, sep='\t', index=False)
        else:
            print("No valid data to write to TSV.")
    else:
        print("No valid data to write to TSV.")

    with open(counter_file, 'w') as cf:
        cf.write(str(record_counter))

# Function to process the low PID log
def main():
    input_dir = 'Flu_fasta'
    output_dir = 'Flu_fasta/output'
    log_file = 'blast_log.txt'
    low_pid_file = 'low_pid_log.txt'
    blast_result_file = 'blast_results.txt'
    key_file_path = 'id_key.tsv'
    counter_file = 'record_counter.txt'
    output_fasta_path = 'filtered_combined_sequences.fasta'
    output_tsv_path = 'filtered_combined_output.tsv'

    os.makedirs(output_dir, exist_ok=True)

    # Initialize the log files
    with open(log_file, 'w') as lf, open(low_pid_file, 'w') as lpf, open(blast_result_file, 'w') as brf:
        lf.write(f'{datetime.datetime.now()} - Main script started\n')
        lpf.write(f'{datetime.datetime.now()} - Starting new run\n')
        brf.write(f'{datetime.datetime.now()} - BLAST results\n')

    all_errors = []
    all_results = []
    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith('.fasta'):
            with open(log_file, 'a') as lf:
                lf.write(f'{datetime.datetime.now()} - Processing file: {fasta_file}\n')
            errors, results = process_fasta_file(
                os.path.join(input_dir, fasta_file), low_pid_file, log_file, blast_result_file
            )
            all_errors.extend(errors)
            all_results.extend(results)
            with open(log_file, 'a') as lf:
                lf.write(f'{datetime.datetime.now()} - Finished processing file: {fasta_file}\n')

    with open(log_file, 'a') as lf:
        for error in all_errors:
            lf.write(error)
        lf.write(f'{datetime.datetime.now()} - Processing BLAST results\n')

    compile_and_filter_sequences(
        input_dir, low_pid_file, output_fasta_path, output_tsv_path, key_file_path, counter_file
    )

    with open(log_file, 'a') as lf:
        lf.write(f'{datetime.datetime.now()} - Main script finished\n')

if __name__ == "__main__":
    main()

