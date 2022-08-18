import os
import sys
import datetime
import threading
import rich_click as click
from collections import defaultdict

from runner import qd_start

from tools.helpers import (
    setup_logger,
    get_config,
    make_samplesheet,
from tools.slims import (
    SlimsSample,
    slims_records_from_sec_analysis,
    find_runtag_from_fastqs,
    find_fastq_paths,
    find_attached_bioinfo_objects,
    update_bioinformatics_record)


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
                name=sample_id,
            )
        )

    # Start all samples in parallel
    finished_samples = []
    for t in threads:
        logger.info(f"{t.name.split('_')[0]} - Starting the runner")
        t.start()
    for u in threads:  # Waits for all threads to finish
        u.join()
        logger.info(f"{u.name.split('_')[0]} - Completed the runner")
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
    all_rnaseq_samples = slims_records_from_sec_analysis(186)
    # 29 = WOPR
    # 186 = QD-RNA

    # Find all DNA objects marked for QD-RNAseq pipeline as secondary analysis
    # They also need to either be missing any bioinformatics objects set as QD-rnaseq,
    # or have a bioinformatics object set as QD-rnaseq with the secAnalysisState set to 'novel'
    rnaseq_samples = defaultdict(dict)
    for sample in all_rnaseq_samples:
        # Create the slimssample object
        slimsinfo = SlimsSample(sample)
        if not slimsinfo.fastqs:  # If no fastqs, skip
            continue

        for fastq in slimsinfo.fastqs:
            # Skip fastq objects set to not include or no bioinformatics
            if fastq.cntn_cstm_doNotInclude.value == True:
                continue
            elif fastq.cntn_cstm_noBioinformaticsObjects.value == True:
                continue

            # Find all bioinformatics objects for this fastq object with correct secondary analysis
            #rnaseq_samples[sample] = {}
            try:
                attached_bioinfo_object = find_attached_bioinfo_objects(slimsinfo.bioinformatics, fastq.pk(), 186)
            except Exception as e:
                logger.error(e)

            # Take care of case where there is no bioinformatics object
            if attached_bioinfo_object == False:
                # Create a bionformatics object and set it to "running"
                bioinfo_fields = {
                    'cntn_cstm_secondaryAnalysis': [186],
                    'cntn_cstm_SecondaryAnalysisState': 'running',
                }
                # Create a bioinformatics object and save it
                new_bioinfo_record = slimsinfo.add_bioinformatics(fastq.pk(), fields=bioinfo_fields)
                rnaseq_samples[sample]['bioinformatics'] = new_bioinfo_record

            # Set existing bioinfo record to 'running' and keep track of it
            elif attached_bioinfo_object.cntn_cstm_SecondaryAnalysisState.value == 'novel':
                new_bioinfo_record = update_bioinformatics_record(attached_bioinfo_object, fields={'cntn_cstm_SecondaryAnalysisState': 'running'})
                rnaseq_samples[sample]['bioinformatics'] = new_bioinfo_record

            else:  # State != novel
                continue

            # Save all the relevant fastq paths for this sample
            if len(rnaseq_samples[sample]) > 0:
                all_fastq_paths = find_fastq_paths(slimsinfo.fastqs)
                rnaseq_samples[sample]['fastq'] = all_fastq_paths

    if len(rnaseq_samples) > 0:
        logger.info(f"Found {len(selected_rnaseq_samples.keys())} sample(s) marked for QD-RNAseq pipeline. Starting the wrapper.")
    else:
        sys.exit(0)

    ### --- Loop over each record --- ###
    runner_samples = {} # Store all samplesheet paths per sample in a dict for later use
    for sample in rnaseq_samples:
        ### --- Collect information about the sample --- ###
        # Create a new SlimsSample object
        logger.info(f"{sample} - Extracting SLIMS information.")
        slimsinfo = SlimsSample(sample)

        # Get the run tag and combine with name
        runtag = find_runtag_from_fastqs(slimsinfo.fastqs)
        sample_id = f"{sample}_{runtag}"

        ### --- Generate a samplesheet from the information gathered --- ###
        outdir = os.path.join(config.get("general", "output_dir"), sample_id)
        logger.info(f"{sample} - Generating samplesheet.csv in {outdir}.")
        os.makedirs(outdir, exist_ok=True)

        strandedness = config.get("general", "strandedness")

        sample_ss_path = make_samplesheet(sample_id, slimsinfo.fastqs, strandedness, outdir)
        runner_samples[sample_id] = sample_ss_path

    ### --- Start a runner for each sample --- ###
    completed_samples = start_runner_threads(runner_samples, logger)

    ### --- Add analysed samples to previously analysed samples file --- ###
    for sample in completed_samples:
        with open(config.get("general", "previously_analysed"), "a") as f:
            f.write(sample + "\n")

    ### --- Completed wrapper --- ###
    logger.info(f"Completed the QD-RNAseq wrapper for {len(completed_samples)} sample(s).")

if __name__ == "__main__":
    main()
