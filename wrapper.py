import os
import sys
import click
import datetime
from tools.helpers import (
    setup_logger,
    get_config,
    dir_to_samplesheet,
    sanitize_fastqdir,
)


@click.command()
@click.option(
    "--fastqdir", required=True, help="Path to input directory of fastq files"
)
@click.option("--outdir", required=True, help="Path to output directory")
@click.option("--strandedness", default="reverse", help="Strandedness of seq libraries")
def main(fastqdir, outdir, strandedness):
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


if __name__ == "__main__":
    main()
