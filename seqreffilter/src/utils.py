import pysam

def parse_bam_by_reads(bam_file):
    """
    parse a bam file, and yield all hit for a read one time
    The bam file should sorted by reads name (samtools sort -n), yield one read one time
    """
    bf = pysam.AlignmentFile(bam_file, 'r')
    read_tmp_record = []
    for r in bf.fetch(until_eof=True):
        if len(read_tmp_record) == 0 or r.query_name == read_tmp_record[-1].query_name:
            read_tmp_record.append(r)
            continue
        else:
            yield read_tmp_record
            read_tmp_record = [r]
    if len(read_tmp_record) != 0:
        yield read_tmp_record
    bf.close()