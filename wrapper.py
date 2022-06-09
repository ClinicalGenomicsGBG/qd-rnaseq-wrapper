import os
import sys
import click
import datetime
import threading
import subprocess
from tools.helpers import (
    setup_logger,
    get_config,
    dir_to_samplesheet,
    sanitize_fastqdir,
    build_rnaseq_command,
    build_rnafusion_command,
)


@click.command()
@click.option(
    "--fastqdir", required=True, help="Path to input directory of fastq files"
)
@click.option("--outdir", required=True, help="Path to output directory")
@click.option("--strandedness", default="reverse", help="Strandedness of seq libraries")
@click.option("--testrun", is_flag=True, help="Run with nf-core test data")
@click.option("--skip-rnaseq", is_flag=True, help="Skip the nf-core/rnaseq pipeline")
@click.option("--skip-rnafusion", is_flag=True, help="Skip the nf-core/rnafusion pipeline")
def main(fastqdir, outdir, strandedness, testrun, skip_rnaseq, skip_rnafusion):
    # Set up the logger function
    now = datetime.datetime.now()
    logfile = os.path.join(outdir, "QD-rnaseq" + now.strftime("%y%m%d_%H%M%S") + ".log")
    # logger = setup_logger("qd-rnaseq", logfile)
    logger = setup_logger("qd-rnaseq")  ## TODO, remove to enable logpath
    logger.info("Starting the RNAseq pipepline wrapper.")

    # Read in the config
    logger.info(f"Reading parameters from config.ini")
    config = get_config()

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

    # Start pipelines
    # Set up the threading
    def call_script(args):
        subprocess.call(args)
    threads = []

    # Build the rnaseq command and add to threads
    if not skip_rnaseq:
        logger.info("Starting the nf-core/rnaseq pipeline")
        rnaseq_command = build_rnaseq_command(config, outdir, ss_path, testrun)
        threads.append(threading.Thread(target=call_script, args=rnaseq_command))

    # Build the rnafusion command and add to threads
    if not skip_rnafusion:
        logger.info("Starting the nf-core/rnafusion pipeline")
        rnafusion_command = build_rnafusion_command(config, outdir, ss_path, testrun)
        threads.append(threading.Thread(target=call_script, args=rnafusion_command))

    # Start both pipelines in parallel
    for t in threads:
        t.start()
    for u in threads: # Waits for all threads to finish
        u.join()


if __name__ == "__main__":
    main()
