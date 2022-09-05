import os
import sys
import rich_click as click
import datetime
from tools.helpers import (
    setup_logger,
    get_config,
    dir_to_samplesheet,
    verify_fastqdir,
    report_results,
)
from tools.processing import (
    build_rnaseq_command,
    build_rnafusion_command,
    start_pipe_processes,
)
from tools.slims import (
    update_bioinformatics_record,
)

@click.group()
def cli():
    pass

@cli.command()
@click.option(
    "--fastqdir",
    help="Path to input directory of fastq files",
)
@click.option(
    "--outdir",
    help="Path to output directory",
)
@click.option(
    "--sample-name",
    help="Name of sample to be used in output. Defaults to fastqdir basename",
)
@click.option(
    "--ss_path",
    help="Path to input samplesheet to use for analysis",
)
@click.option(
    "--strandedness",
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
@click.option(
    "--skip-report",
    is_flag=True,
    help="Skips moving output data to reporting folder",
)
def cli_qd_start(
    fastqdir,
    outdir,
    sample_name,
    ss_path,
    strandedness,
    testrun,
    skip_rnaseq,
    skip_rnafusion,
    save_reference,
    skip_report,
):
    qd_start(
        fastqdir=fastqdir,
        outdir=outdir,
        sample_name=sample_name,
        ss_path=ss_path,
        strandedness=strandedness,
        testrun=testrun,
        skip_rnaseq=skip_rnaseq,
        skip_rnafusion=skip_rnafusion,
        save_reference=save_reference,
        skip_report=skip_report)

def qd_start(
    fastqdir: str = None,
    outdir: str = None,
    sample_name: str = None,
    ss_path: str = None,
    strandedness: str = None,
    testrun: bool = False,
    skip_rnaseq: bool = False,
    skip_rnafusion: bool = False,
    save_reference: bool = False,
    skip_report: bool = False,
    logger = None,
) -> None:

    # Read in the config
    config = get_config()

    # Set up the logger function
    if not logger:
        now = datetime.datetime.now()
        logdir = config.get("general", "wrapper_log_dir") # TODO, this should not be a log on wrapper level
        os.makedirs(logdir, exist_ok=True)
        logfile = os.path.join(
            logdir,
            "QD-rnaseq-wrapper_" + now.strftime("%y%m%d_%H%M%S") + ".log",
        )
        logger = setup_logger("qd-rnaseq", logfile)

    # If testrun, skip fastq handling
    if not testrun:
        #Sanity check
        if fastqdir == None and ss_path == None:
            logger.error("Please provide either a fastqdir or a samplesheet")
            sys.exit(1)

        # Make samplesheet from fastq dir
        if not ss_path:
            logger.info(f"No samplesheet provided. Creating samplesheet.csv from {fastqdir}")
            scriptpath = config.get("general", "fastq_to_ss_path")
            try:
                verify_fastqdir(fastqdir)
                ss_path = dir_to_samplesheet(scriptpath, fastqdir, strandedness)
            except Exception as e:
                logger.error(e)
                sys.exit(1)

    else:
        sample_name = 'testrun'
        ss_path = ''  # No samplesheet in testruns
        logger.info("Running pipelines with test data")

    # Set the samplename and output directory. Default to fastqdir basename
    if not sample_name:
        sample_name = os.path.basename(os.path.normpath(fastqdir))
    if outdir is None:
        outdir = os.path.join(config.get("general", "output_dir"), sample_name)
        os.makedirs(outdir, exist_ok=True)


    # Empty list for storing which pipes to start
    pipe_commands = {}

    # Build the rnaseq command and add to processes
    if not skip_rnaseq:
        rnaseq_command = build_rnaseq_command(
            config,
            outdir,
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
            ss_path,
            testrun,
        )
        # Merge the rnafusion command with the other commands
        pipe_commands = pipe_commands | rnafusion_command

    # Start the pipelines in separate threads
    finished_pipes = start_pipe_processes(sample_name, pipe_commands, logger)

    # Move selected files to report dir
    if skip_report:
        logger.info("Reporting disabled, skipping copy to report dir")
    else:
        report_dir = os.path.join(config.get("general", "report_dir"), sample_name)
        transferred_files = report_results(finished_pipes, outdir, sample_name, config)
        for pipe, transfers in transferred_files.items():
            logger.info(f"Moved {transfers} files generated by {pipe}")

    #TODO, set runstate in the bioinformatics object
    #This is just psueod-code, not yet tested
    #if completed prooerly:
        # update_bioinformatics_record(bioinformatics_object, fields={'cntn_cstm_SecondaryAnalysisState': 'novel'})
    #elif completed_incorrec:
        # update_bioinformatics_record(bioinformatics_object, fields={'cntn_cstm_SecondaryAnalysisState': 'error'})


    # Pipeline completed
    logger.info("Completed the RNAseq wrapper workflow")

if __name__ == "__main__":
    cli_qd_start()
