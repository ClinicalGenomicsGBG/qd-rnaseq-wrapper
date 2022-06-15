import os
import logging
from configparser import ConfigParser
import subprocess


def setup_logger(name, log_path=None):
    """Enables the logging to be setup correctly"""
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


def dir_to_samplesheet(scriptpath: str, fastqdir: str, strandedness: str) -> str:
    """
    Executes the nf-core fastq_dir_to_samplesheet.py script on a given dir with fastq files
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


def sanitize_fastqdir(fastqdir: str):
    """
    Checks all files in a given fastq dir that assumptions about naming is met.
    For now this only checks file endings, but could be expanded to check for samples with more than 2 files.
    """
    for filename in os.listdir(fastqdir):
        if not filename.endswith(".fastq.gz"):
            if filename.lower() == "samplesheet.csv":
                continue
            raise Exception(
                f"Found a file not ending with .fastq.gz in fastq dir: {filename}"
            )


def build_rnaseq_command(
    config, outdir: str, logdir: str, ss_path: str, testrun=False, save_reference=False
) -> list:
    """
    Create a list with the nextflow commands to send to subprocess

    :param config: configparser object with configurations
    :param outdir: Path to output directory
    :param logdir: Path to directory to store logs
    :param ss_path: Path to samplesheet.csv
    :param: testrun: Set to true to run test data
    :param: save_reference: Set to save the downloaded reference genomes in output dir
    :return: List of list with all components of the command
    """
    rnaseq_command = ["nextflow"]

    # Place of logs
    rnaseq_command.append("-log")
    rnaseq_command.append(os.path.join(logdir, "rnaseq", "rnaseq.log"))

    # Path to main.nf
    rnaseq_command.append("run")
    rnaseq_command.append(config.get("nextflow", "rnaseq"))

    # Execution report
    rnaseq_command.append("-with-report")
    rnaseq_command.append(os.path.join(logdir, "rnaseq", "rnaseq-execution.html"))

    # If custom config, use it
    if config.get("nextflow", "custom_config"):
        rnaseq_command.append("-c")
        rnaseq_command.append(config.get("nextflow", "custom_config"))

    # Profile to use
    rnaseq_command.append("-profile")
    if testrun:
        rnaseq_command.append(config.get("nextflow", "test_profile"))
    else:
        rnaseq_command.append(config.get("nextflow", "profile"))


    # Pre-downloaded references, or download from iGenomes
    if config.has_section("rnaseq-references") and not save_reference:
        for option in config.options("rnaseq-references"):
            rnaseq_command.append(f"--{option}")
            rnaseq_command.append(config.get("rnaseq-references", option))
    else:
        rnaseq_command.append("--genome")
        rnaseq_command.append(config.get("nextflow", "genome"))

    # Input samplesheet.csv
    if not testrun:
        rnaseq_command.append("--input")
        rnaseq_command.append(ss_path)

    # Outdir
    rnaseq_command.append("--outdir")
    rnaseq_command.append(os.path.join(outdir, "rnaseq"))

    if save_reference:
        rnaseq_command.append("--save_reference")

    return rnaseq_command

def build_rnafusion_command(config, outdir: str, logdir: str, ss_path: str, testrun=False) -> list:
    """
    Create a list with the nextflow commands to send to subprocess

    :param config: configparser object with configurations
    :param outdir: Path to output directory
    :param logdir: Path to directory to store logs
    :param ss_path: Path to samplesheet.csv
    :param: testrun: Set to true to run test data
    :return: List of list with all components of the command
    """
    rnafusion_command = ["nextflow"]

    # Place of logs
    rnafusion_command.append("-log")
    rnafusion_command.append(os.path.join(logdir, "rnafusion", "rnafusion.log"))

    # Path to main.nf
    rnafusion_command.append("run")
    rnafusion_command.append(config.get("nextflow", "rnafusion"))

    # Execution report
    rnafusion_command.append("-with-report")
    rnafusion_command.append(
        os.path.join(logdir, "rnafusion", "rnafusion-execution.html")
    )

    # If custom config, use it
    if config.get("nextflow", "custom_config"):
        rnafusion_command.append("-c")
        rnafusion_command.append(config.get("nextflow", "custom_config"))

    # Profile to use
    rnafusion_command.append("-profile")
    if testrun:
        rnafusion_command.append(config.get("nextflow", "test_profile"))
    else:
        rnafusion_command.append(config.get("nextflow", "profile"))

    # Run alla tools
    rnafusion_command.append("--all")

    # Path to dependencies
    rnafusion_command.append("--genomes_base")
    rnafusion_command.append(config.get("nextflow", "dependencies_fusion"))

    # Use filtered fusionreport fusions for fusioninspector
    rnafusion_command.append("--fusioninspector_filter")

    # Input samplesheet.csv
    if not testrun:
        rnafusion_command.append("--input")
        rnafusion_command.append(ss_path)

    # Outdir
    rnafusion_command.append("--outdir")
    rnafusion_command.append(os.path.join(outdir, "rnafusion"))

    return rnafusion_command
