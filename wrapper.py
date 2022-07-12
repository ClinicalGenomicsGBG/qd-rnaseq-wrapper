import rich_click as click
from tools.slims import SlimsSample, translate_slims_info, samples_from_sec_analysis, fastq_paths

@click.command()
@click.option(
    "--log",
    help="Path to wrapper log file. Default value set in config.ini",
)
def main():
    # Find all slims records marked for QD-RNAseq pipeline as secondary analysis
    rnaseq_samples = samples_from_sec_analysis(186)
    # 29 = WOPR
    # 186 = QD-RNA

    # Loop over each record and gather relevant information
    # Start the runner in a separate thread for each sample
    for sample, record in rnaseq_samples.items():
        # Create a new SlimsSample object
        slimsinfo = SlimsSample(sample)


        # Save fastq paths for forward and reverse in separate lists
        fastqs = fastq_paths(slimsinfo.fastqs)
        fastq_forwards = []
        fastq_reverses = []
        for fastq_run in fastqs:
            fastq_forwards.append(fastq_run[1][0])
            fastq_reverses.append(fastq_run[1][1])

        #print(f'Forwards: {" ".join(fastq_forwards)}')
        #print(f'Reverses: {" ".join(fastq_forwards)}')

        # Translate the information from the SLIMS database into a dictionary
        info = translate_slims_info(record)


if __name__ == "__main__":
    main()