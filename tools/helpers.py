import os
import sys
import csv
import glob
import shutil
import logging
from configparser import ConfigParser
import subprocess
import threading
from .slims import fastq_paths


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


def build_rnaseq_command(
    config,
    outdir: str,
    ss_path: str,
    testrun=False,
    save_reference=False,
) -> dict[str, list[str]]:
    """
    Create a list with the nextflow commands to send to subprocess

    :param config: configparser object with configurations
    :param outdir: Path to output directory
    :param ss_path: Path to samplesheet.csv
    :param testrun: Set to true to run test data
    :param save_reference: Save the downloaded reference genomes in output dir
    :return: Dict of list with all components of the command
    """
    rnaseq_command = ["nextflow"]

    # Place of logs
    rnaseq_command.append("-log")
    rnaseq_command.append(os.path.join(outdir, "logs", "rnaseq.log"))

    # Path to main.nf
    rnaseq_command.append("run")
    rnaseq_command.append(config.get("rnaseq", "main"))

    # Execution report
    rnaseq_command.append("-with-report")
    rnaseq_command.append(
        os.path.join(
            outdir,
            "logs",
            "rnaseq-execution.html",
        )
    )

    # If custom config, use it
    if config.has_option("nextflow", "custom_config"):
        rnaseq_command.append("-c")
        rnaseq_command.append(config.get("nextflow", "custom_config"))
    elif os.path.isfile(
        os.path.join(
            os.path.dirname(sys.modules["__main__"].__file__),
            "qd-wrapper_nextflow.config",
        )
    ):
        rnaseq_command.append("-c")
        rnaseq_command.append(
            os.path.join(
                os.path.dirname(sys.modules["__main__"].__file__),
                "qd-wrapper_nextflow.config",
            )
        )

    # Profile to use
    rnaseq_command.append("-profile")
    if testrun:
        rnaseq_command.append(
            config.get(
                "nextflow",
                "test_profile",
            )
        )
    else:
        rnaseq_command.append(
            config.get(
                "nextflow",
                "profile",
            )
        )

    # Outdir
    rnaseq_command.append("--outdir")
    rnaseq_command.append(os.path.join(outdir, "rnaseq"))

    # If testrun, all input is now completed
    if testrun:
        return {"nf-core/rnaseq": rnaseq_command}

    # Pre-downloaded references, or download from iGenomes
    if config.has_section("rnaseq-references") and not save_reference:
        for option in config.options("rnaseq-references"):
            #Deal with which index to use separately
            if option == "star_index" or option == "rsem_index":
                continue
            rnaseq_command.append(f"--{option}")
            rnaseq_command.append(config.get("rnaseq-references", option))
    else:
        rnaseq_command.append("--genome")
        rnaseq_command.append(config.get("rnaseq", "genome"))

    # Aligner to use
    rnaseq_command.append("--aligner")
    rnaseq_command.append(config.get("rnaseq", "aligner"))
    # Add salmon if not using the star_salmon option
    if config.get("rnaseq", "aligner") != "star_salmon":
        rnaseq_command.append("--pseudo_aligner")
        rnaseq_command.append("salmon")

    # Set the index to use
    aligner = config.get("rnaseq", "aligner")

    if aligner == "star_salmon":
        rnaseq_command.append("--star_index")
        rnaseq_command.append(config.get("rnaseq-references", "star_index"))
    elif aligner == "star_rsem":
        rnaseq_command.append("--rsem_index")
        rnaseq_command.append(config.get("rnaseq-references", "rsem_index"))
    else:
        raise Exception(f"Aligner {aligner} not supported.")

    # Input samplesheet.csv
    rnaseq_command.append("--input")
    rnaseq_command.append(ss_path)

    if save_reference:
        rnaseq_command.append("--save_reference")

    return {"nf-core/rnaseq": rnaseq_command}


def build_rnafusion_command(
    config,
    outdir: str,
    ss_path: str,
    testrun=False,
) -> dict[str, list[str]]:
    """
    Create a list with the nextflow commands to send to subprocess

    :param config: configparser object with configurations
    :param outdir: Path to output directory
    :param ss_path: Path to samplesheet.csv
    :param testrun: Set to true to run test data
    :return: Dict of list with all components of the command
    """
    rnafusion_command = ["nextflow"]

    # Place of logs
    rnafusion_command.append("-log")
    rnafusion_command.append(os.path.join(outdir, "logs", "rnafusion.log"))

    # Path to main.nf
    rnafusion_command.append("run")
    rnafusion_command.append(config.get("rnafusion", "main"))

    # Execution report
    rnafusion_command.append("-with-report")
    rnafusion_command.append(os.path.join(outdir, "logs", "rnafusion-execution.html"))

    # If custom config, use it
    if config.has_option("nextflow", "custom_config"):
        rnafusion_command.append("-c")
        rnafusion_command.append(config.get("nextflow", "custom_config"))
    elif os.path.isfile(
        os.path.join(
            os.path.dirname(sys.modules["__main__"].__file__),
            "qd-wrapper_nextflow.config",
        )
    ):
        rnafusion_command.append("-c")
        rnafusion_command.append(
            os.path.join(
                os.path.dirname(sys.modules["__main__"].__file__),
                "qd-wrapper_nextflow.config",
            )
        )

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
    rnafusion_command.append(config.get("rnafusion", "dependencies_fusion"))

    # Use filtered fusionreport fusions for fusioninspector
    rnafusion_command.append("--fusioninspector_filter")

    # Set fusionreport tool cutoff
    rnafusion_command.append("--fusionreport-tool-cutoff")
    rnafusion_command.append(config.get("rnafusion", "fusionreport_tool_cutoff"))

    # Set the readlength to use
    rnafusion_command.append("--read_length")
    rnafusion_command.append(config.get("rnafusion", "readlength"))

    # Input samplesheet.csv
    if not testrun:
        rnafusion_command.append("--input")
        rnafusion_command.append(ss_path)

    # Outdir
    rnafusion_command.append("--outdir")
    rnafusion_command.append(os.path.join(outdir, "rnafusion"))

    return {"nf-core/rnafusion": rnafusion_command}


def start_pipe_threads(sample_name:str, pipe_dict: dict, logger) -> list:
    """
    Takes a dict where keys are pipeline names and values are a
    list containing all parts of a command, and starts them in a
    separate thread using subprocess.call

    :param sample_name: Name of the sample, for logging purposes
    :param pipe_dict: Dict of lists where each one is a command to run
    :param logger: Logger object to write logs to
    """
    # Set up the threading
    def call_script(args):
        subprocess.call(args)

    # Create the threading object
    threads = []
    for pipe, command in pipe_dict.items():
        threads.append(
            threading.Thread(
                target=call_script,
                args=[command],
                name=pipe,
            )
        )

    # Start both pipelines in parallel
    finished_pipes = []
    for t in threads:
        logger.info(f"{sample_name.split('_')[0]} - Starting the {t.name} pipeline")
        t.start()
    for u in threads:  # Waits for all threads to finish
        u.join()
        logger.info(f"{sample_name.split('_')[0]} - Completed the {u.name} pipeline")
        finished_pipes.append(u.name)

    return finished_pipes


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

def make_samplesheet(sample, fastqs, strandedness: str, outdir: str) -> str:
    """
    Takes a list of fastq files and a strandedness and creates a samplesheet.csv
    file with the correct information.

    :param sample: Name of the sample
    :param fastqs: Object with fastq file paths
    :param strandedness: Strandedness of the library
    :param outfile: Path to the output file
    :return: Path to the samplesheet.csv file
    """
    # Get all fastq paths
    fastqs = fastq_paths(fastqs)

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