import os
import sys
import uuid
import shutil
import datetime
import threading
import rich_click as click

from runner import qd_start

from tools.helpers import (
    setup_logger,
    get_config,
    make_samplesheet)
from tools.slims import (
    SlimsSample,
    translate_slims_info,
    samples_from_sec_analysis,
    find_runtag_from_fastqs)

def start_runner_threads(sample_dict: dict, logger) -> list:
    """
    Takes a dict containing sample names and samlesheet paths and starts an instance
    of the runner in a separate threads.

    :param sample_dict: Dict of samples with their samplesheet paths
    :param logger: Logger object to write logs to
    :return: List of finished samples
    """
    threads = []
    for sample_id, ss_path in sample_dict.items():
        qd_start_kwargs = {'sample_name': sample_id, 'ss_path': ss_path, 'logger': logger}

        threads.append(
            threading.Thread(
                target=qd_start,
                kwargs=qd_start_kwargs,
                name=sample_id.split('_')[0],
            )
        )

    # Start all samples in parallel
    finished_samples = []
    for t in threads:
        logger.info(f"{t.name} - Starting the runner")
        t.start()
    for u in threads:  # Waits for all threads to finish
        u.join()
        logger.info(f"{u.name}- Completed the runner")
        finished_samples.append(u.name)

    return finished_samples

@click.command()
@click.option(
    "--logdir",
    help="Path to directory of wrapper log file. Default value set in config.ini",
)
@click.option(
    "--cleanup",
    help="Set to remove any files in /tmp on completion of wrapper",
    is_flag=True,
)
def main(logdir: str, cleanup: bool):
    ### --- Read in the config --- ###
    config = get_config()

    #### --- Set up the logger function --- ###
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
    rnaseq_samples = samples_from_sec_analysis(186)
    # 29 = WOPR
    # 186 = QD-RNA
    if len(rnaseq_samples) > 0:
        logger.info(f"Found {len(rnaseq_samples)} sample(s) marked for QD-RNAseq pipeline. Starting the wrapper.")
    else:
        sys.exit(0)

    ### --- Make a temp dir for current invocation of the wrapper --- ###
    random_uuid = str(uuid.uuid4().hex)
    temp_dir = os.path.join(config.get("general", "temp_dir"), random_uuid)
    temp_dir = '/tmp/qd-testing' #TODO, remove this line
    try:
        os.makedirs(temp_dir)
        logger.info(f"Creating a temp dir for files in {temp_dir}")
    except Exception:
        logger.error(f"Could not create temp directory: {temp_dir}")
        sys.exit(1)

    ### --- Loop over each record --- ###
    runner_samples = {} # Store all samplesheet paths per sample in a dict for later use
    for sample, record in rnaseq_samples.items():
        ### --- Collect information about the sample --- ###
        # Create a new SlimsSample object
        logger.info(f"{sample} - Extracting SLIMS information.")
        slimsinfo = SlimsSample(sample)

        # Translate the information from the SLIMS database into a dictionary
        info = translate_slims_info(record)

        # Get the runtag for the sample, combine with sample name
        try:
            runtag = find_runtag_from_fastqs(slimsinfo.fastqs)
            logger.info(f"{sample} - Setting {runtag} as run tag.")
            sample_id = f"{sample}_{runtag}"
        except Exception as e:
            logger.error(e)
            sys.exit(1)
            # TODO, this should mark this sample as failed, and send an email. Not break the whole thing

        ### --- Generate a samplesheet from the information gathered --- ###
        outdir = os.path.join(config.get("general", "output_dir"), sample_id)
        logger.info(f"{sample} - Generating samplesheet.csv in {outdir}.")
        os.makedirs(outdir, exist_ok=True)

        strandedness = config.get("general", "strandedness")

        sample_ss_path = make_samplesheet(sample_id, slimsinfo.fastqs, strandedness, outdir)
        runner_samples[sample_id] = sample_ss_path

    ### --- Start a runner for each sample --- ###
    completed_samples = start_runner_threads(runner_samples, logger)

    ### --- Clean up temporary files --- ###
    if cleanup:
        logger.info(f"Cleaning up temporary files in {temp_dir}")
        shutil.rmtree(temp_dir)

    ### --- Completed wrapper --- ###
    logger.info(f"Completed the QD-RNAseq wrapper for {len(completed_samples)} sample(s).")

if __name__ == "__main__":
    main()
