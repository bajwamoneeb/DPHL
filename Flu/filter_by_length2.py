import os
import datetime
from Bio import SeqIO
import pandas as pd

# Define the length cutoffs based on GISAID's criteria
LENGTH_CUTOFFS = {
    "HA_H1": (1701, 1785),
    "HA_H2": (1701, 1785),  
    "HA_H3": (1701, 1785),
    "NA_N1": (1408, 1467),
    "NA_N2": (1408, 1467),
    "NA_N3": (1408, 1467),  
    "MP": (800, 1027),
    "NS": (700, 980),
    "PA": (2152, 2300),
    "PB1": (2275, 2341),
    "PB2": (2265, 2334),
    "NP": (1474, 1569),
    "B_HA": (1700, 1850),
    "B_MP": (950, 1150),
    "B_NA": (1350, 1500),
    "B_NS": (900, 1100),
    "B_PA": (2100, 2300),
    "B_PB1": (2300, 2500),
    "B_PB2": (2300, 2500),
    "B_NP": (1500, 1700)
}

# Function to extract the correct segment name
def extract_segment_name(record_id):
    parts = record_id.split("_")
    if len(parts) >= 3 and parts[1] == 'B':
        return "B_" + parts[2]
    elif len(parts) >= 3 and parts[1] == 'A':
        if parts[2] in ["HA", "NA"]:
            return parts[2] + "_" + parts[3] if len(parts) > 3 else parts[2]
        return parts[2]
    else:
        return parts[-1]

# Function to merge duplicate columns based on base column names
def merge_duplicate_columns(df):
    columns_to_merge = {}
    
    for col in df.columns:
        base_col = col.split('.')[0]  # Remove .1, .2, etc.
        if base_col in columns_to_merge:
            columns_to_merge[base_col].append(col)
        else:
            columns_to_merge[base_col] = [col]
    
    for base_col, col_list in columns_to_merge.items():
        if len(col_list) > 1:
            df[base_col] = df[col_list].apply(lambda row: row.dropna().iloc[0] if not row.dropna().empty else None, axis=1)
            df.drop(columns=col_list[1:], inplace=True)
        else:
            df.rename(columns={col_list[0]: base_col}, inplace=True)
    
    return df

# Function to filter sequences based on length and track logs
def filter_fasta_by_length(input_fasta, output_fasta, length_cutoffs, log_file, id_key_file, sample_status):
    id_key_df = pd.read_csv(id_key_file, sep='\t')  # Read the ID key file

    with open(input_fasta) as in_handle, open(output_fasta, "w") as out_handle, open(log_file, "w") as log_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            segment = extract_segment_name(record.id)
            original_id = id_key_df.loc[id_key_df['New_ID'] == record.id, 'Original_ID'].values[0]
            file_prefix = original_id.split("-")[0]

            if file_prefix not in sample_status:
                sample_status[file_prefix] = "PASS"

            if segment in length_cutoffs:
                min_length, max_length = length_cutoffs[segment]
                if min_length <= len(record.seq) <= max_length:
                    out_handle.write(f">{record.id}\n{str(record.seq)}\n")
                else:
                    log_handle.write(f"{file_prefix}\tRemoved {record.id} - Length: {len(record.seq)} (cutoff: {min_length}-{max_length})\n")
                    sample_status[file_prefix] = "FAIL"
            else:
                log_handle.write(f"Segment {record.id} does not have a specified length cutoff.\n")
                sample_status[file_prefix] = "FAIL"

# Function to update the PID log file with the File_Prefix and combine logs
def update_pid_log_file(pid_log_file, updated_pid_log_file, sample_status):
    if not os.path.exists(pid_log_file):
        print(f"PID log file {pid_log_file} not found, skipping PID log update.")
        return

    pid_logs = {}
    with open(pid_log_file) as pid_file:
        for line in pid_file:
            # Skip the starting log line if it exists
            if "Starting new run" in line or "No PIDs were below the cutoff" in line:
                continue
            if "Segment" in line:
                parts = line.strip().split(" ")
                file_prefix = parts[6]
                log = " ".join(parts[3:])
                if file_prefix in pid_logs:
                    pid_logs[file_prefix].append(log)
                else:
                    pid_logs[file_prefix] = [log]
                sample_status[file_prefix] = "FAIL"

    with open(updated_pid_log_file, "w") as updated_pid_file:
        updated_pid_file.write("File_Prefix\tPID_log\n")
        for file_prefix, logs in pid_logs.items():
            combined_logs = "; ".join(logs)
            updated_pid_file.write(f"{file_prefix}\t{combined_logs}\n")

def combine_length_logs(length_log_file, combined_length_log_file):
    length_logs = {}
    with open(length_log_file) as length_file:
        for line in length_file:
            # Skip empty or invalid lines
            if "\t" not in line:
                continue
            try:
                file_prefix, log = line.strip().split("\t", 1)
                if file_prefix in length_logs:
                    length_logs[file_prefix].append(log)
                else:
                    length_logs[file_prefix] = [log]
            except ValueError:
                print(f"Skipping malformed line: {line.strip()}")

    with open(combined_length_log_file, "w") as combined_length_file:
        combined_length_file.write("File_Prefix\tLength_filter_log\n")
        for file_prefix, logs in length_logs.items():
            combined_logs = "; ".join(logs)
            combined_length_file.write(f"{file_prefix}\t{combined_logs}\n")

# Function to combine the PID and length filter logs
def combine_logs(updated_pid_log_file, combined_length_log_file, combined_log_file, sample_status):
    if not os.path.exists(updated_pid_log_file):
        pid_log_df = pd.DataFrame(columns=["File_Prefix", "PID_log"])
    else:
        pid_log_df = pd.read_csv(updated_pid_log_file, sep='\t')

    length_log_df = pd.read_csv(combined_length_log_file, sep='\t')

    combined_df = pd.merge(pid_log_df, length_log_df, on='File_Prefix', how='outer')
    combined_df.rename(columns={'File_Prefix': 'entity:flu_id'}, inplace=True)
    
    combined_df['HA_NA_QC'] = combined_df.apply(lambda row: 'FAIL' if isinstance(row['Length_filter_log'], str) and ('NA' in row['Length_filter_log'] or 'HA' in row['Length_filter_log']) else 'PASS', axis=1)
    
    combined_df['HA_NA_QC'] = combined_df['HA_NA_QC'].fillna('PASS')

    for sample in sample_status:
        if sample_status[sample] == "PASS" and sample not in combined_df['entity:flu_id'].values:
            combined_df = pd.concat([combined_df, pd.DataFrame({'entity:flu_id': [sample], 'HA_NA_QC': ['PASS']})], ignore_index=True)
        elif sample_status[sample] == "FAIL" and sample not in combined_df['entity:flu_id'].values:
            combined_df = pd.concat([combined_df, pd.DataFrame({'entity:flu_id': [sample], 'HA_NA_QC': ['FAIL']})], ignore_index=True)

    combined_df.to_csv(combined_log_file, sep='\t', index=False)

# Function to update the TSV file based on the filtered FASTA file
# Also logs samples missing or failing HA or NA length check
def update_tsv_file(fasta_file, input_tsv, output_tsv, sample_status, length_log_file):
    # Get the filtered IDs from the FASTA file
    filtered_ids = {record.id for record in SeqIO.parse(fasta_file, "fasta")}
    
    # Read the TSV file
    df = pd.read_csv(input_tsv, sep='\t')

    # Merge duplicate columns if any
    df = merge_duplicate_columns(df)

    # Open the log file to append entries
    with open(length_log_file, "a") as log_handle:
        # Check for missing or invalid HA and NA segments
        for index, row in df.iterrows():
            file_prefix = row['File_Prefix']
            ha_valid = row['Seq_Id (HA)'] in filtered_ids if pd.notnull(row['Seq_Id (HA)']) else False
            na_valid = row['Seq_Id (NA)'] in filtered_ids if pd.notnull(row['Seq_Id (NA)']) else False

            # If either HA or NA is missing or invalid, mark the sample as FAIL
            if not ha_valid or not na_valid:
                sample_status[file_prefix] = "FAIL"
                log_handle.write(f"{file_prefix}\tRemoved sample - Missing or failed HA or NA segment\n")
                df.at[index, 'Seq_Id (HA)'] = None
                df.at[index, 'Seq_Id (NA)'] = None
        
        # Drop rows where HA or NA segments are missing or invalid
        rows_to_remove = df[df[['Seq_Id (HA)', 'Seq_Id (NA)']].isnull().any(axis=1)].index
        df.drop(rows_to_remove, inplace=True)

    # Write the updated TSV file
    df.to_csv(output_tsv, sep='\t', index=False)

def main():
    input_fasta_file = 'filtered_combined_sequences.fasta'
    input_tsv_file = 'filtered_combined_output.tsv'
    id_key_file = 'id_key.tsv'
    pid_log_file = 'low_pid_log.txt'
    output_fasta_file = 'filtered_combined_sequences_2.fasta'
    output_tsv_file = 'filtered_tsv_length_cutoff.tsv'
    length_log_file = 'length_filter_script_log.txt'
    combined_length_log_file = 'length_log.txt'
    updated_pid_log_file = 'updated_pid_log.txt'
    combined_log_file = 'combined_log.txt'

    sample_status = {}

    # Perform the filtering of the FASTA file
    filter_fasta_by_length(input_fasta_file, output_fasta_file, LENGTH_CUTOFFS, length_log_file, id_key_file, sample_status)

    # Update the TSV file with the new filtering and failure conditions
    update_tsv_file(output_fasta_file, input_tsv_file, output_tsv_file, sample_status, length_log_file)

    # Update the PID log file, accommodating potential empty or irrelevant content
    update_pid_log_file(pid_log_file, updated_pid_log_file, sample_status)

    # Combine the length and PID logs
    combine_length_logs(length_log_file, combined_length_log_file)

    for prefix in sample_status.keys():
        if sample_status[prefix] == "PASS":
            print(f"Sample {prefix} passed all filters.")

    combine_logs(updated_pid_log_file, combined_length_log_file, combined_log_file, sample_status)

if __name__ == "__main__":
    main()
