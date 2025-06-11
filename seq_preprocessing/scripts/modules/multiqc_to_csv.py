from sys import argv
import pandas as pd
import os

# This script extracts a specific column from a tsv (meant for the MultiQC table output)
# And adds the values as a new column in an existing csv table, while matching the values in the sample column
# This is used to collect info on read length and read count from MultiQC and append it on separate files after every step
# Data_table must have a column "Sample" that matches a column "sample" in out_table
def multiqc_to_csv(data_table, out_table, data_col, out_col):
  # Read data tables
  df_data = pd.read_csv(data_table, sep='\t')
  df_out = pd.read_csv(out_table)
   
  # Some checks
  # Check if data_col exists in df_data, if not, break
  if data_col not in df_data.columns:
    print(f"Column '{data_col}' does not exist in {os.path.basename(data_table)}")
    return
  # Then check if out_col exists in df_out, if it exists already it will be overwritten
  elif out_col in df_out.columns:
    print(f"Column '{out_col}' already exists in {os.path.basename(out_table)}. Overwriting...")
    
    # Delete column
    df_out = df_out.drop(out_col, axis=1)
    
  else:
    print(f"Adding data to column '{out_col}'")
  
  # Merge data based on 'Sample' column
  merged_df = pd.merge(df_out, df_data[['Sample', data_col]], left_on='sample', right_on='Sample', how='outer')
  merged_df['sample'] = merged_df['sample'].fillna(merged_df['Sample'])
  merged_df = merged_df.drop("Sample", axis=1)
  
  # Rename the new column
  merged_df = merged_df.rename(columns={data_col: out_col})
  merged_df.to_csv(out_table, index=False)

if __name__ == "__main__":
  input_table = argv[1]
  output_table = argv[2]
  input_col = argv[3]
  output_col = argv[4]
  multiqc_to_csv(input_table, output_table, input_col, output_col)
