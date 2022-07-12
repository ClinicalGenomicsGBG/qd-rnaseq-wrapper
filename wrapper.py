import os
import datetime
import rich_click as click

from tools.helpers import (
    setup_logger,
    get_config,)
from tools.slims import (
    SlimsSample,
    translate_slims_info,
    samples_from_sec_analysis,
    fastq_paths,)

@click.command()
@click.option(
    "--logdir",
    help="Path to directory of wrapper log file. Default value set in config.ini",
)
def main(logdir):
    ### --- Read in the config --- ###
    config = get_config()

    #### ------- Set up the logger function --- ###
    now = datetime.datetime.now()
    if not logdir:
        logdir = config.get("general", "wrapper_log_dir")

    os.makedirs(logdir, exist_ok=True)

    logfile = os.path.join(
        logdir,
        "QD-rnaseq-wrapper_" + now.strftime("%y%m%d_%H%M%S") + ".log",
        )
    logger = setup_logger("qd-rnaseq", logfile)
    logger.info("Starting the RNAseq pipepline wrapper.")

    ### --- Find all slims records marked for QD-RNAseq pipeline as secondary analysis --- ###
    rnaseq_samples = samples_from_sec_analysis(186)
    # 29 = WOPR
    # 186 = QD-RNA

    ### --- Loop over each record --- ###
    # Gather relevant information for each sample
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