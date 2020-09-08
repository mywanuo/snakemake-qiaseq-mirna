import pandas as pd

def merge_read_counts(input_files, output_file):
    """
    Merges sample-specific RNA counts into one .tsv file.
    First column contains the RNA identifier and rest of
    the columns represent individual samples.
    """

    # Loop through each input file.
    for i in range(0, len(input_files)):

        # Print status message.
        print("Reading counts from {}...".format(input_files[i]))

        # Import count table as pd.DataFrame,
        # use first column as an index.
        df = pd.read_table(input_files[i], index_col="Molecule")

        # Extract sample ID from file path.
        # Rename second column according to sample ID.
        sample_id = input_files[i].split("/")[1]
        df = df.rename(columns={"Count":sample_id})

        if i == 0:
            # If this is the first count table to iterate over, use it as a
            # template for merged count table.
            merged_df = df
        else:
            # If iterated count table is not the first one, join table with
            # the pre-existing merged table.
            merged_df = merged_df.join(df, how="outer")

    # Fill empty values (no expression) with zero.
    merged_df = merged_df.fillna(0)

    # Write merged pd.DataFrame to tab-delimited file.
    merged_df.to_csv(output_file[0], sep="\t")

# Call method.
merge_read_counts(snakemake.input, snakemake.output)