import os
import sys
from multiprocessing import Process
import subprocess

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
    processes = []
    for pipe, command in pipe_dict.items():
        threads.append(
            Process(
                target=call_script,
                args=[command],
                name=pipe,
            )
        )

    # Start both pipelines in parallel
    finished_pipes = []
    for p in processes:
        logger.info(f"{sample_name} - Starting the {p.name} pipeline")
        p.start()
    for p in processes:  # Waits for all threads to finish
        p.join()
        logger.info(f"{sample_name} - Completed the {p.name} pipeline")
        finished_pipes.append(p.name)

    return finished_pipes