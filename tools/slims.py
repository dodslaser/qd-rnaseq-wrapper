import os
import json
from datetime import datetime

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction
from slims.content import Status

import passwords

# Initial set up of the SLIMS instance
instance = 'qdRNAseq_query'
url = passwords.slims.url
user = passwords.slims.user
password = passwords.slims.password
slims = Slims(instance, url, user, password)


class SlimsSample:
    """
    Class containing the sample data.
    """
    def __init__(self, content_id):
        self.content_id = content_id
        self._dna = None
        self._fastqs = []
        self._bioinformatics = []

    @property
    def dna(self):
        if not self._dna:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.content_id))
                                  .add(equals('cntn_fk_contentType', 6)))

            if len(records) > 1:
                raise Exception('More than 1 DNA somehow.')

            if records:
                self._dna = records[0]

        return self._dna


    @property
    def fastqs(self):
        if not self._fastqs:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.content_id))
                                  .add(equals('cntn_fk_contentType', 22)))
            self._fastqs = records

        return self._fastqs


    @property
    def bioinformatics(self):
        if not self._bioinformatics:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.content_id))
                                  .add(equals('cntn_fk_contentType', 23)))
            self._bioinformatics = records

        return self._bioinformatics

    def add_bioinformatics(self, original_content_pk, fields={}):
        if not self.fastqs:
            raise Exception('Can not add bioinformatics without a parent fastq record.')

        fields['cntn_id'] = self.content_id
        fields['cntn_fk_contentType'] = 23
        fields['cntn_status'] = Status.PENDING.value
        fields['cntn_fk_location'] = 83
        fields['cntn_fk_originalContent'] = original_content_pk
        fields['cntn_fk_user'] = ''

        self._bioinformatics.append(slims.add('Content', fields))

        return self._bioinformatics

    def refresh(self):
        self._dna = None
        self._fastq = None
        self._bioinformatics = None


def translate_slims_info(record) -> dict:
    """
    Takes a SLIMS record and translates relevant data into more easily handled dict

    :param record: SLIMS record
    :return: Dict containing relevant data about the sample
    """
    sample_name = record.cntn_id.value

    pipeline_pks = record.cntn_cstm_secondaryAnalysis.value
    pcr = record.cntn_cstm_pcr.value

    investigator = 'CGG'  # NOTE: Needed?

    department_record = slims.fetch_by_pk('ReferenceDataRecord', record.cntn_cstm_department.value)
    department = department_record.rdrc_name.value  # Format KK
    responder_records = [slims.fetch_by_pk('ReferenceDataRecord', pk) for
                         pk in department_record.rdrc_cstm_responder.value]
    responder_emails = [rec.rdrc_cstm_email.value for rec in responder_records]

    is_research = record.cntn_cstm_research.value
    research_project = record.cntn_cstm_researchProject.value

    is_priority = True if record.cntn_cstm_priority.value else False

    gender = record.gender.value

    is_trio = record.cntn_cstm_trio.value
    trio_id = record.cntn_cstm_trioID.value
    trio_role = record.cntn_cstm_trioRole.value

    tertiary_analysis = record.cntn_cstm_tertiaryAnalysis.value

    master = {
        'content_id': sample_name,
        'investigator': investigator,
        'department': department,
        'responder_mails': responder_emails,
        'is_research': is_research,
        'research_project': research_project,
        'gender': gender,
        'is_priority': is_priority,
        'pcr': pcr,
        'is_trio': is_trio,
        'trio_id': trio_id,
        'trio_role': trio_role,
        'secondary_analysis': pipeline_pks,
        'tertiary_analysis': tertiary_analysis
    }
    return master


def slims_records_from_sec_analysis(primary_key: int, content_type: int = 6) -> list:
    """
    Fetch all objects from slims given a secondary analysis key and content type

    :param slims: A slims connection object
    :param primary_key: Primary key integer of secondary analysis
    :param content_type Key of content type, default is 6 (DNA)
    :return: List with all records marked with the secondary analysis key
    """

    records = slims.fetch("Content",conjunction()
                          .add(equals("cntn_cstm_secondaryAnalysis", primary_key))
                          .add(equals("cntn_fk_contentType", content_type)))

    return records


def find_fastq_paths(fastq_records) -> list:
    """
    Takes a slims object with fastq content and returns a list of paths to the fastq files
    with the number of reads per pair

    :param fastq_records: slims object with fastq content
    :return: list of fastq objects with the number of reads and a tuple of the paths to the fastq files
    """
    reads_fastq_pairs = []  # NOTE: All reads and fastq_paths tuple pairs
    for fastq_object in fastq_records:
        # Skip if set as not include
        if fastq_object.cntn_cstm_doNotInclude.value == True:
            continue

        demuxer_info_json = json.loads(fastq_object.cntn_cstm_demuxerSampleResult.value)
        if not demuxer_info_json:
            raise MissingDemuxerInfoError(f'A fastq object without demuxer info was found: {sample_name}')

        fastq_paths = demuxer_info_json['fastq_paths']
        reads = demuxer_info_json['total_reads']

        for fastq_path in fastq_paths:
            if not os.path.exists(fastq_path):
                raise MissingFastqError(f'A fastq path from Slims does not exist locally: {fastq_path}')

        reads_fastq_pairs.append((reads, fastq_paths))  # Tuple
    return reads_fastq_pairs

def find_runtag_from_fastqs(fastq_objects) -> str:
    """
    Takes a slims object with fastq content and returns the run tag. This is sort of hacky,
    and it would be better if this information was directly available in the SlimsObject

    :param fastq_objects: slims object with fastq content
    :return: Run tag (YYMMDD_flowcallID)
    """
    # First get all runtags associated with all fastq files and store in list
    runtag_list = []
    for fastq_object in fastq_objects:
        # Grab objects run_tag
        remote_run_tag = fastq_object.cntn_cstm_runTag.value
        runtag_list.append(remote_run_tag)

    # Loop over the list, and figure out which is the most recent runtag
    runtag_newest = None # Initialise a newest object to compare against

    for runtag in runtag_list:
        runtag_date = runtag.split('_')[0]
        runtag_date_dt = datetime.strptime(runtag_date, "%y%m%d")

        if runtag_newest == None:
            runtag_newest = runtag
        else:
            runtag_newest_dt = datetime.strptime(runtag_newest.split('_')[0], "%y%m%d")
            if runtag_date_dt > runtag_newest_dt:
                runtag_newest = runtag

    if runtag_newest == None:
        raise Exception('No run tag found in fastq objects')
    else:
        return runtag_newest


def find_secanalysis_and_state(bioinformatics_record, sample_name) -> tuple[int, str]:
    """
    Takes a slims object with bioinformatics content and returns a tuple containing
    the secondary analysis ID and secondary analysis state

    :param bioinformatics_records: slims bioinformatics object
    :param sample_name: Sample name (for logging purposes)
    :return: Tuple with secondary analysis ID and secondary analysis state
    """
    # Note, this is a bit hacky since right now a bioinfo object could have multiple secondary analyses.
    # This means that cntn_cstm_secondaryAnalysis is a list
    sec_analysis = bioinformatics_record.cntn_cstm_secondaryAnalysis.value[0]
    sec_analysis_state = bioinformatics_record.cntn_cstm_SecondaryAnalysisState.value
    if not sec_analysis:
        raise Exception(f'A bioinformatics object without secondary analysis info was found: {sample_name}')
    elif not sec_analysis_state:
        raise Exception(
            f'A bioinformatics object without secondary analysis state was found: {sample_name}')

    return sec_analysis, sec_analysis_state


def fetch_bioinformatics_record(content_id: str, secondary_analysis: int) -> 'Records':
    """
    Fetch a bioinformatics record from slims given a content ID, content type, and secondary analysis ID

    :param content_id: Content ID
    :param secondary_analysis: Secondary analysis ID
    :return: SLIMS bioinformatics object
    """
    records = slims.fetch("Content", conjunction()
                          .add(equals("cntn_id", content_id))
                          .add(equals("cntn_fk_contentType", 23))  # 23 equals bioinformatics
                          .add(equals("cntn_cstm_secondaryAnalysis", secondary_analysis)))

    if len(records) > 1:
        raise IOError(f'Multiple bioinformatics objects found for {content_id}')
    else:
        return records


def fetch_fastq_records(content_id: str,) -> 'Record':
    """
    Fetch fastq records from slims given a content ID

    :param content_id: Content ID
    :return: SLIMS fastq object
    """
    records = slims.fetch("Content", conjunction()
                          .add(equals("cntn_id", content_id))
                          .add(equals("cntn_fk_contentType", 22)))  # 22 equals fastq

    return records


def find_derived_bioinfo_objects(bioinformatics_records, original_content_pk: int, secondaryAnalysis: int) -> 'Record':
    """
    Find all bioinformatics objects attached to a given originalContent and secondaryAnalysis

    :param bioinformatics_records: SLIMS bioinformatics records
    :param original_content_pk: Original content primary key
    :param secondaryAnalysis: Secondary analysis key
    :return: Bioinformatics objects
    """

    attached_bioinfo_objects = []
    for bioinfo_obj in bioinformatics_records:

        if (
                bioinfo_obj.cntn_fk_originalContent.value == original_content_pk
                and secondaryAnalysis in bioinfo_obj.cntn_cstm_secondaryAnalysis.value
        ):
            attached_bioinfo_objects.append(bioinfo_obj)

    if len(attached_bioinfo_objects) > 1:
        raise Exception(f"Content {original_content_pk} has more than one bioinformatics object attached with secondary analysis {secondaryAnalysis}")

    if len(attached_bioinfo_objects) == 0:
        return False
    else:
        return attached_bioinfo_objects[0]


def update_bioinformatics_record(bioinformatics_record, fields={}) -> 'Record':
    """
    Update a bioinformatics record in slims given a dictionary of fields to update

    :param bioinformatics_records: SLIMS bioinformatics records
    :param fields: Dictionary of fields to update
    :return: SLIMS bioinformatics object
    """

    updated_record = bioinformatics_record.update(fields)

    return updated_record
