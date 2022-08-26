import os
import glob
import shutil
import logging
from configparser import ConfigParser
import threading
import subprocess


def setup_logger(name, log_path=None):
    """Enables the logging to be setup correctly"""
    # Remove any existing logging, from SLIMS package for instance.
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    stream_handle = logging.StreamHandler()
    stream_handle.setLevel(logging.DEBUG)
    stream_handle.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(stream_handle)

    if log_path:
        file_handle = logging.FileHandler(log_path, "a")
        file_handle.setLevel(logging.DEBUG)
        file_handle.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(file_handle)

    return logger


def get_config():
    """
    Read in a config.ini file and return the object. Assumes it's located in
    the wrapper root folder and called config.ini
    """
    converters = {"list": lambda x: [i.strip() for i in x.split(",")]}
    config = ConfigParser(converters=converters)

    root_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    config_path = os.path.join(root_path, "config.ini")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"No config found at {config_path}")

    config.read(config_path)
    return config


def dir_to_samplesheet(
    scriptpath: str,
    fastqdir: str,
    strandedness: str,
) -> str:
    """
    Executes the nf-core fastq_dir_to_samplesheet.py script on a
    given dir with fastq files
    :param scriptpath: Path to fastq_dir_to_samplesheet.py
    :param fastqdir: Path to directory with fastq files
    :param strandedness: Strandedness of library
    :return: Path to the generated samplesheet.csv
    """
    ss_path = os.path.join(fastqdir, "samplesheet.csv")
    command_args = [
        "python",
        scriptpath,
        fastqdir,
        ss_path,
        "--strandedness",
        strandedness,
    ]

    try:
        subprocess.run(command_args)
    except:
        raise Exception(f"Could not execute {scriptpath}.")

    return ss_path


def verify_fastqdir(fastqdir: str) -> None:
    """
    Checks all files in a given fastq dir that assumptions about naming is met.
    For now this only checks file endings, but could be expanded to check for
    samples with more than 2 files.
    """
    for filename in os.listdir(fastqdir):
        if not filename.endswith(".fastq.gz"):
            if filename.lower() == "samplesheet.csv":
                continue
            raise Exception(
                f"Found a file not ending with .fastq.gz in fastq dir: {filename}"
            )


def report_results(finished_pipes: list, outdir: str, sample_name: str, config) -> dict:
    """
    Takes a list of finished pipelines and checks which output files from this
    pipeline should be included in the final report (from config file). Copies these to the
    specified report directory in the config. Returns a dict with the
    number of copied files per pipeline

    :param finished_pipes: List of finished pipelines
    :param outdir: Path from where to move files (pipeline outdir)
    :param sample_name: Name of the sample
    :param config: Configparser object with configurations
    :return: Dict with the number of copied files per pipeline
    """
    # Read in where stuff should be copied to
    report_dir = os.path.join(config.get("general", "report_dir"), sample_name)
    igv_dir = os.path.join(config.get("general", "igv_dir"), sample_name)

    # Setup a dict for returning how many files
    copied_files = {}

    # Copy files for each pipeline used
    for pipeline in finished_pipes:
        # Remove chars before the slash ("nf-core/")
        pipeline = pipeline.split("/")[-1]

        # Read in the config file for this pipeline
        for option in config.options(f"report-{pipeline}"):

            # Unpack all options
            for file_refex in config.get(f"report-{pipeline}", option).split(','):

                # Special case for stringtie as this is nested under the name of the aligner
                if option == 'stringtie':
                    aligner = config.get('rnaseq', 'aligner')
                    search_path = os.path.join(outdir, pipeline, aligner, option, file_refex)
                else:
                    search_path = os.path.join(outdir, pipeline, option, file_refex)

                # for all files in path
                for file in glob.glob(search_path):
                    # Copy the file to the correct location
                    destination = os.path.join(
                            report_dir,
                            option,
                            os.path.basename(file))
                    os.makedirs(os.path.dirname(destination), exist_ok=True)

                    if file.endswith(('.bam','.bai','.bigWig')): #IGV files
                        destination = os.path.join(
                            igv_dir,
                            os.path.basename(file))
                        os.makedirs(os.path.dirname(destination),
                                    exist_ok=True)

                        shutil.copy(file, destination)

                    elif os.path.isdir(file):
                        # Copy the directory
                        shutil.copytree(file, destination, dirs_exist_ok=True)
                    else:
                        # Copy the file
                        shutil.copy(file, destination)

                    # Add to the dict
                    if pipeline in copied_files:
                        copied_files[pipeline] += 1
                    else:
                        copied_files[pipeline] = 1

    return copied_files

def make_samplesheet(sample: str, fastqs: list, strandedness: str, outdir: str) -> str:
    """
    Takes a list of fastq files and a strandedness and creates a samplesheet.csv
    file with the correct information.

    :param sample: Name of the sample
    :param fastqs: List with fastq file paths
    :param strandedness: Strandedness of the library
    :param outfile: Path to the output file
    :return: Path to the samplesheet.csv file
    """

    # Path to samplesheet.csv
    ss_path = os.path.join(outdir, "samplesheet.csv")

    # Create the samplesheet
    with open(ss_path, 'w') as f:
        f.write("sample,fastq_1,fastq_2,strandedness\n")
        for fastq in fastqs:
            f.write(f"{sample},{fastq[1][0]},{fastq[1][1]},{strandedness}\n")

    return ss_path

def read_previous_samples_file(config) -> list:
    """
    From the config, read the file containing the previous samples and return a
    list of all samples.

    :param config: Configparser object with configurations
    :return: List of sample names
    """
    previous_samples_file = config.get("general", "previously_analysed")

    with open(previous_samples_file, "r") as prev:
        previous_samples = [line.rstrip() for line in prev]

    return previous_samples