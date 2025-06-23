from sys import argv
import pandas as pd
import os
import re

# This script takes as input the idxstats output from samtools, which shows coverage per contig
# and a table linking contigs to reference genomes, and generates a summarised coverage table per genome (instead of per contig)
def summarise_mapping_stats(sample_name, idxstats, contig_table_path):
    idxstats.columns=["contig", "length", "mapped", "unmapped"] # Add colnames
    sum_mapped=idxstats["mapped"].sum() # Get sum of mapped sequences
    # Add ref name to idxstats table
    merged_df=pd.merge(idxstats, contig_table, on="contig")
    # Summarise mapped reads per group
    summary_df=merged_df.groupby('reference').agg({'mapped': 'sum'}).T
    summary_df.insert(0, "sample", sample_name)
    # Get host species name
    names=summary_df.columns.tolist()
    filtered_names = [name for name in names if name not in ["sample", "Homo sapiens", "PhiX"]]
    if filtered_names: # Replace host species name with "host"
        host_species = filtered_names[0]
        summary_df.rename(columns={host_species: 'host'}, inplace=True)
    else: # If there is no host name, because sample is a blank or control, just return NA under "host"
        host_species = "NA"  # Or any other value you prefer
        summary_df["host"]="NA"
    print("Total mapped reads for " + sample_name + ": " + str(sum_mapped))
    # Replace it with host
    return(summary_df)

def loop_idxstats(suffix, contig_table):
    # List idxstats files
    txt_files = [file for file in os.listdir() if file.endswith(suffix)]
    # Initiate an empty table to populate
    output_table=pd.DataFrame(columns=["sample", "host", "Homo sapiens", "PhiX"])
    for idxstats_path in txt_files:
        # Print sample name
        sample_name=idxstats_path.split('_mapped_idxstats.txt')[0]
        # Read idxstats table
        idxstats = pd.read_csv(idxstats_path, sep='\t', header=None)
        # Run summarise function
        df=summarise_mapping_stats(sample_name, idxstats, contig_table)
        # Print some stats
        formatted_values = []
        output_table=pd.concat([output_table, df])
        row = df[df["sample"] == sample_name]
        for column in output_table.columns[1:]:
            formatted_value = f"{column}: {row[column].values[0]}"
            formatted_values.append(formatted_value)
        print(", ".join(formatted_values))
    return(output_table)

if __name__ == "__main__":
    # Get args from command line
    suffix = argv[1]
    contig_table_path = argv[2]
    # Load contig table
    contig_table = pd.read_csv(contig_table_path, sep=",").drop_duplicates()
    # Remove ">" from the beginning of contigs
    contig_table["contig"]=contig_table["contig"].apply(lambda x: x.split(">")[1])
    # Run loop function
    output_table=loop_idxstats(suffix, contig_table)
    # get outpath
    outdir=os.getcwd()
    output_table.to_csv(os.path.join(outdir, "genome_mapping_stats.tsv"), sep='\t', index=False) 
