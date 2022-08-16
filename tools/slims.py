import os
import json
from datetime import datetime

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction
from slims.content import Status

import passwords

# Initial set up of the SLIMS instance
instance = 'wopr_query'
url = passwords.slims.url
user = passwords.slims.user
password = passwords.slims.password
slims = Slims(instance, url, user, password)


class SlimsSample:
    """
    Class containing the sample data.
    """
    def __init__(self, sample_name, run_tag=''):
        self.sample_name = sample_name
        self.run_tag = run_tag

        self._dna = None
        self._fastq = None
        self._fastqs = []
        self._bioinformatics = None

    @property
    def dna(self):
        if not self._dna:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 6)))

            if len(records) > 1:
                raise Exception('More than 1 DNA somehow.')

            if records:
                self._dna = records[0]

        return self._dna

    @property
    def fastq(self):
        if not self.run_tag:
            raise Exception('Can not fetch fastq without a set run tag.')

        if not self._fastq:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22))
                                  .add(equals('cntn_cstm_runTag', self.run_tag)))

            if len(records) > 1:
                raise Exception('More than 1 fastq somehow.')

            if records:
                self._fastq = records[0]

        return self._fastq

    @property
    def fastqs(self):
        if not self._fastqs:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22)))
            self._fastqs = records

        return self._fastqs

    @property
    def bioinformatics(self):
        if not self.run_tag:
            raise Exception('Can not fetch fastq without a set run tag.')

        if not self._bioinformatics:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 23))
                                  .add(equals('cntn_cstm_runTag', self.run_tag)))

            if len(records) > 1:
                raise Exception('More than 1 bioinformatics somehow.')

            if records:
                self._bioinformatics = records[0]

        return self._bioinformatics

    def add_bioinformatics(self, fields={}):
        if self.bioinformatics:
            raise Exception('Bioinformatics content already exists. Can not add.')
        if not self.fastq:
            raise Exception('Can not add bioinformatics without a parent fastq record.')
        if not self.run_tag:
            raise Exception('Can not add bioinformatics without a set run tag.')

        fields['cntn_id'] = self.sample_name
        fields['cntn_fk_contentType'] = 23
        fields['cntn_cstm_runTag'] = self.run_tag
        fields['cntn_status'] = Status.PENDING.value
        fields['cntn_fk_location'] = 83
        fields['cntn_fk_originalContent'] = self.fastq.cntn_pk.value
        fields['cntn_fk_user'] = ''

        self._bioinformatics = slims.add('Content', fields)

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

def samples_from_sec_analysis(primary_key: int) -> dict:
    """
    Fetch all records from slims given a secondary analysis tag

    :param primary_key: Primary key integer of secondary analysis
    :return: Dictionary of sample names and their corresponding slims records
    """

    records = slims.fetch("Content", equals("cntn_cstm_secondaryAnalysis", primary_key))

    samples = {}

    for record in records:
        sample_name = record.column('cntn_id').value
        samples[sample_name] = record

    return samples


def find_fastq_paths(fastq_records) -> list:
    """
    Takes a slims object with fastq content and returns a list of paths to the fastq files
    with the number of reads per pair

    :param fastq_records: slims object with fastq content
    :return: list of fastq objects with the number of reads and a tuple of the paths to the fastq files
    """
    reads_fastq_pairs = []  # NOTE: All reads and fastq_paths tuple pairs
    for fastq_object in fastq_records:
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
