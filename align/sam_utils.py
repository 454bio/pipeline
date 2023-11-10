import pysam

header = { 'HD': {'VN': '1.0'},
           'SQ': [{'LN': 5386, 'SN': 'phiX174'}]}

class SamUtils:
    def __init__(self, filename:str = 'out.sam'):
        self.fp = pysam.AlignmentFile(filename, "w", header=header)

    def __del__(self):
        self.fp.close()

    def AddRead(self, read:str, pos:int, ref:str, read_name:str = 'read'):
        l = len(read)
        qual = '<'*l # for now just provide default qual
        a = pysam.AlignedSegment()
        a.query_name = read_name
        a.query_sequence = read
        a.reference_id = 0 # index in the SQ header array
        a.reference_start = pos
        a.mapping_quality = 20 # just default qual for now
        a.cigarstring = self.GenCigar(read, ref)
        a.template_length = 5386 # need to allow genome and chrom/ref inputs to include setup in header, etc
        a.query_qualities = pysam.qualitystring_to_array(qual)

        self.fp.write(a)

    def GenCigar(self, read, ref):
        # very basic, we will just look for mismatches and mark the positions
        cigar = ''
        for i in range(len(read)):
            if read[i] != ref[i]:
                cigar += str(i+1) + 'M'
        return cigar

