import os
import sys
import rich_click as click
import datetime
from tools.helpers import (
    setup_logger,
    get_config,
    dir_to_samplesheet,
    sanitize_fastqdir,
    build_rnaseq_command,
    build_rnafusion_command,
    start_pipe_threads,
)


@click.command()
@click.option(
    "--fastqdir",
    help="Path to input directory of fastq files",
)
@click.option(
    "--outdir",
    required=True,
    help="Path to output directory",
)
@click.option(
    "--strandedness",
    default="reverse",
    help="Strandedness of seq libraries",
)
@click.option(
    "--testrun",
    is_flag=True,
    help="Run with nf-core test data",
)
@click.option(
    "--skip-rnaseq",
    is_flag=True,
    help="Skip the nf-core/rnaseq pipeline",
)
@click.option(
    "--skip-rnafusion",
    is_flag=True,
    help="Skip the nf-core/rnafusion pipeline",
)
@click.option(
    "--save-reference",
    is_flag=True,
    help="Save the genome references in outfolder",
)
def main(
    fastqdir: str,
    outdir: str,
    strandedness: str,
    testrun: bool,
    skip_rnaseq: bool,
    skip_rnafusion: bool,
    save_reference: bool,
) -> None:
    # Set up the logger function
    now = datetime.datetime.now()
    logdir = os.path.join(outdir, "logs")
    os.makedirs(logdir, exist_ok=True)
    logfile = os.path.join(
        logdir,
        "QD-rnaseq-wrapper_" + now.strftime("%y%m%d_%H%M%S") + ".log",
    )
    logger = setup_logger("qd-rnaseq", logfile)
    logger.info("Starting the RNAseq pipepline wrapper.")

    # Read in the config
    logger.info(f"Reading parameters from config.ini")
    config = get_config()

    # If testrun, skip fastq handling
    if not testrun:
        # Make sure fastqdir looks good
        try:
            sanitize_fastqdir(fastqdir)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

        # Make samplesheet from fastq dir
        logger.info(f"Creating samplesheet.csv from {fastqdir}")
        scriptpath = config.get("general", "fastq_to_ss_path")
        try:
            ss_path = dir_to_samplesheet(scriptpath, fastqdir, strandedness)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
    else:
        logger.info("Running pipelines with test data")

    # Set the samplename and output directory. Default to fastqdir basename
    if testrun:
        sample_name = 'testrun'
        ss_path = '' # No samplesheet in testruns
    elif not sample_name:
        sample_name = os.path.basename(os.path.normpath(fastqdir))
    if outdir is None:
        outdir = os.path.join(config.get("general", "output_dir"), sample_name)

    os.makedirs(outdir, exist_ok=True)
    logger.info(f"Output directory: {outdir}")

    # Get the name of the sample
    # TODO, this should be read from the samplesheet
    # and then looped over to run the pipelines
    logger.info(f"Processing sample: {sample_name}")

    # Empty list for storing which pipes to start
    pipe_commands = {}

    # Build the rnaseq command and add to threads
    if not skip_rnaseq:
        rnaseq_command = build_rnaseq_command(
            config,
            outdir,
            logdir,
            ss_path,
            testrun,
            save_reference,
        )
        # Merge the rnaseq command with the other commands
        pipe_commands = pipe_commands | rnaseq_command

    # Build the rnafusion command and add to threads
    if not skip_rnafusion:
        rnafusion_command = build_rnafusion_command(
            config,
            outdir,
            logdir,
            ss_path,
            testrun,
        )
        # Merge the rnafusion command with the other commands
        pipe_commands = pipe_commands | rnafusion_command

    # Start the pipelines in separate threads
    start_pipe_threads(pipe_commands, logger)

    # Pipeline completed
    logger.info("Completed the RNAseq wrapper workflow")


if __name__ == "__main__":
    main()
