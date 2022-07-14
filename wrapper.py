import os
import datetime
import rich_click as click

from tools.helpers import (
    setup_logger,
    get_config,
    make_samplesheet,)
from tools.slims import (
    SlimsSample,
    translate_slims_info,
    samples_from_sec_analysis,
    find_runtag_from_fastqs)

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


    ### --- Find all slims records marked for QD-RNAseq pipeline as secondary analysis --- ###
    logger.info("Starting the RNAseq pipeline wrapper.")
    rnaseq_samples = samples_from_sec_analysis(186)
    # 29 = WOPR
    # 186 = QD-RNA

    ### --- Loop over each record --- ###
    for sample, record in rnaseq_samples.items():
        ### --- Collect information about the sample --- ###
        # Create a new SlimsSample object
        logger.info(f"Extracting SLIMS information for {sample}.")
        slimsinfo = SlimsSample(sample)

        # Translate the information from the SLIMS database into a dictionary
        info = translate_slims_info(record)

        # Get the runtag for the sample
        try:
            runtag = find_runtag_from_fastqs(slimsinfo.fastqs)
            logger.info(f"Setting {runtag} as run tag for sample {sample}.")
        except Exception as e:
            logger.error(e)
            sys.exit(1)
            # TODO, this should mark this sample as failed, and send an email. Not break the whole thing

        ### --- Generate a samplesheet from the information gathered --- ###
        outdir = os.path.join(config.get("general", "output_dir"), sample)
        logger.info(f"Generating samplesheet.csv for {sample} in {outdir}.")
        os.makedirs(outdir, exist_ok=True)

        strandedness = config.get("general", "strandedness")

        make_samplesheet(sample, slimsinfo.fastqs, strandedness, outdir)


if __name__ == "__main__":
    main()