import os
import click
import datetime
from tools.helpers import setup_logger, get_config

@click.command()
@click.option("--fastqdir", required=True, help="Path to input directory of fastq files"
)
@click.option("--outdir", required=True, help="Path to output directory"
)
def main(fastqdir, outdir):
    # Set up the logger function
    now = datetime.datetime.now()
    logfile = os.path.join(
        outdir, "QD-rnaseq" + now.strftime("%y%m%d_%H%M%S") + ".log"
    )
    #logger = setup_logger("qd-rnaseq", logfile)
    logger = setup_logger("qd-rnaseq") ## TODO, remove to enable logpath
    logger.info("Starting the RNAseq pipepline wrapper.")

    # Read in the config
    logger.info(f"Reading parameters from config.ini")
    config = get_config()

if __name__ == "__main__":
    main()
