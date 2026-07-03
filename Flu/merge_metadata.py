import pandas as pd

def merge_metadata_and_barcodes(barcode_file, metadata_file):
    # Read the barcode file (assuming it's tab-separated)
    barcode_df = pd.read_csv(barcode_file, sep='\t')
    
    # Read the metadata file (since it's a CSV)
    metadata_df = pd.read_csv(metadata_file)

    # Merge barcodes with metadata on 'Molecular Barcodes' and 'Label Id'
    merged_metadata_df = pd.merge(metadata_df, barcode_df, left_on='Label Id', right_on='Molecular Barcodes', how='inner')
    
    # Strip leading/trailing spaces in 'entity:flu_id' to avoid mismatch issues
    merged_metadata_df['entity:flu_id'] = merged_metadata_df['entity:flu_id'].str.strip()

    # Define constant values
    merged_metadata_df['Location'] = 'United States'
    merged_metadata_df['province'] = 'Delaware'
    merged_metadata_df['Host'] = 'Human'
    merged_metadata_df['Originating_Lab_Id'] = '1345'
    
    # Parse the 'Sampled Date' with the known format: MM/DD/YYYY HH:MM:SS AM/PM
    merged_metadata_df['Collection_Date'] = pd.to_datetime(merged_metadata_df['Sampled Date'], format='%m/%d/%Y %I:%M:%S %p', errors='coerce')
    
    # Remove the timestamp and keep only the date in YYYY-MM-DD format
    merged_metadata_df['Collection_Date'] = merged_metadata_df['Collection_Date'].dt.strftime('%Y-%m-%d')
    
    # Create Collection_Month and Collection_Year
    merged_metadata_df['Collection_Month'] = pd.to_datetime(merged_metadata_df['Collection_Date']).dt.month
    merged_metadata_df['Collection_Year'] = pd.to_datetime(merged_metadata_df['Collection_Date']).dt.year

    return merged_metadata_df, barcode_df

def merge_with_output_tsv(merged_metadata_df, barcode_df, output_tsv_file, unmatched_file_prefixes_file):
    # Read the output TSV file (assuming it's tab-separated)
    output_tsv_df = pd.read_csv(output_tsv_file, sep='\t')
    
    # Strip leading/trailing spaces in 'File_Prefix' to avoid mismatch issues
    output_tsv_df['File_Prefix'] = output_tsv_df['File_Prefix'].str.strip()

    # Perform an outer merge to capture all File_Prefixes, including those that don't match
    merged_df = pd.merge(merged_metadata_df, output_tsv_df, left_on='entity:flu_id', right_on='File_Prefix', how='outer', indicator=True)
    
    # Keep 'Molecular Barcodes' in the final output, even for non-matching File_Prefixes
    merged_df['Molecular Barcodes'] = merged_df['Molecular Barcodes'].fillna('')

    # Extract Isolate Name from the first part of 'Seq_Id (HA)' or any other relevant sequence column
    merged_df['Isolate_Name'] = merged_df['Seq_Id (HA)'].apply(lambda x: x.split('_')[0] if pd.notnull(x) else None)

    # Derive Subtype from `Seq_Id (HA)` and `Seq_Id (NA)`
    def derive_subtype(row):
        ha_id = row['Seq_Id (HA)']
        na_id = row['Seq_Id (NA)']
    
    # Check if `Seq_Id (HA)` is not NaN and contains "B"
        if pd.notnull(ha_id) and "B" in ha_id:
            return "B"  # or "B/Victoria" if lineage information is available
    
        # Otherwise, derive the subtype from HA and NA suffixes if both are available
        ha_suffix = ha_id.split('_')[-1] if pd.notnull(ha_id) else ""
        na_suffix = na_id.split('_')[-1] if pd.notnull(na_id) else ""
    
        if ha_suffix and na_suffix:
            return f"{ha_suffix}{na_suffix}"
    
        return None

    merged_df['Subtype'] = merged_df.apply(derive_subtype, axis=1)

    # Remove duplicates based on 'entity:flu_id' and 'File_Prefix' (or all columns if needed)
    merged_df = merged_df.drop_duplicates(subset=['entity:flu_id', 'File_Prefix'])

    # Separate out the non-matching File_Prefixes
    unmatched_df = merged_df[merged_df['_merge'] == 'right_only']

    # Attempt to match the unmatched 'File_Prefix' with 'Molecular Barcodes' in the barcode_df (fbarcodes file)
    unmatched_df = unmatched_df.merge(barcode_df[['Molecular Barcodes', 'entity:flu_id']], left_on='File_Prefix', right_on='entity:flu_id', how='left', suffixes=('', '_from_fbarcodes'))

    # Save unmatched File_Prefixes and their corresponding Molecular Barcodes (if found) to a file
    unmatched_df[['File_Prefix', 'Molecular Barcodes_from_fbarcodes']].to_csv(unmatched_file_prefixes_file, sep=',', index=False)
    print(f'Unmatched File_Prefixes and corresponding Molecular Barcodes saved to {unmatched_file_prefixes_file}')
    
    # Keep only the rows where there is a match ('both') for the final output
    final_merged_df = merged_df[merged_df['_merge'] == 'both'].drop(columns=['_merge'])

    # Add the necessary columns for GISAID
    additional_columns = [
        'Isolate_Id', 'Segment_Ids', 'Antigen_Character', 'Adamantanes_Resistance_geno', 
        'Oseltamivir_Resistance_geno', 'Zanamivir_Resistance_geno', 'Peramivir_Resistance_geno', 
        'Other_Resistance_geno', 'Adamantanes_Resistance_pheno', 'Oseltamivir_Resistance_pheno', 
        'Zanamivir_Resistance_pheno', 'Peramivir_Resistance_pheno', 'Other_Resistance_pheno', 
        'Host_Age', 'Host_Age_Unit', 'Host_Gender', 'Health_Status', 'Note', 'PMID'
    ]

    for col in additional_columns:
        if col not in final_merged_df.columns:
            final_merged_df[col] = ""

    # Define the columns required for GISAID upload
    desired_columns = [
        'Isolate_Name', 'Isolate_Id', 'Segment_Ids', 'Subtype', 'Lineage', 'Passage_History', 
        'Location', 'province', 'sub_province', 'Location_Additional_info', 'Host', 'Host_Additional_info', 
        'Seq_Id (HA)', 'Seq_Id (NA)', 'Seq_Id (PB1)', 'Seq_Id (PB2)', 'Seq_Id (PA)', 'Seq_Id (MP)', 'Seq_Id (NS)', 
        'Seq_Id (NP)', 'Seq_Id (HE)', 'Seq_Id (P3)', 'Submitting_Sample_Id', 'Authors', 'Originating_Lab_Id', 
        'Originating_Sample_Id', 'Collection_Month', 'Collection_Year', 'Collection_Date', 
        'Antigen_Character', 'Adamantanes_Resistance_geno', 'Oseltamivir_Resistance_geno', 
        'Zanamivir_Resistance_geno', 'Peramivir_Resistance_geno', 'Other_Resistance_geno', 
        'Adamantanes_Resistance_pheno', 'Oseltamivir_Resistance_pheno', 'Zanamivir_Resistance_pheno', 
        'Peramivir_Resistance_pheno', 'Other_Resistance_pheno', 'Host_Age', 'Host_Age_Unit', 'Host_Gender', 
        'Health_Status', 'Note', 'PMID'
    ]

    # Keep only the relevant columns in the final output
    final_merged_df = final_merged_df.reindex(columns=desired_columns)

    # Safely remove 'File_Prefix' if it exists
    if 'File_Prefix' in final_merged_df.columns:
        final_merged_df = final_merged_df.drop(columns=['File_Prefix'])

    return final_merged_df

def main():
    barcode_file = 'Flu_fbarcodes.txt' 
    metadata_file = 'FLU_LIMS_ALL.csv'  
    output_tsv_file = 'filtered_tsv_length_cutoff.tsv'  
    output_merged_file = 'final_merged_flu_metadata_output.csv'  
    unmatched_file_prefixes_file = 'unmatched_file_prefixes.csv'  # File to store non-matching File_Prefixes and barcodes
    
    # Step 1: Merge the metadata and barcodes
    merged_metadata_df, barcode_df = merge_metadata_and_barcodes(barcode_file, metadata_file)
    
    # Step 2: Merge with the output TSV, remove duplicates, and capture unmatched File_Prefixes
    final_merged_df = merge_with_output_tsv(merged_metadata_df, barcode_df, output_tsv_file, unmatched_file_prefixes_file)
    
    # Write the final merged data to a CSV file, keeping only the columns needed for GISAID uploads
    final_merged_df.to_csv(output_merged_file, sep=',', index=False)
    print(f'Final merged data written to {output_merged_file}')

if __name__ == "__main__":
    main()
