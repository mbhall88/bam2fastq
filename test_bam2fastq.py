"""Test cases for the bam2fastq.py script"""
import unittest
import pysam
import bam2fastq
from io import StringIO

HEADER = pysam.AlignmentHeader().from_text("""@HD	VN:1.0	SO:coordinate
@SQ	SN:R00000042	LN:5231428	AS:gi|26111730|gb|AE014075.1|	SP:Ecol
@RG	ID:824f45e8-37f3-4cb9-8a05-63f0b7c9b959	PL:ILLUMINA	PU:160129_D00417_0381_AHJ2VGBCXX_2	LB:VAU2662A45	DT:2016-01-31T00:00:00+0000	SM:H125100459	CN:WTCHG
@RG	ID:7f568ff7-e0f6-4a55-ad17-6fe778ed8f83	PL:ILLUMINA	PU:160129_D00417_0381_AHJ2VGBCXX_1	LB:VAU2662A45	DT:2016-01-31T00:00:00+0000	SM:H125100459	CN:WTCHG
@CO	ID:stampy	TM:Mon, 07 Mar 2016 17:58:14 GMT	WD:/tmp/usecase3938954872417714817dir	HN:gel-pipeline3	UN:compass
@CO	ID:stampy	TM:Mon, 07 Mar 2016 20:24:40 GMT	WD:/tmp/usecase3938954872417714817dir	HN:gel-pipeline3	UN:compass
@CO	PN:stampy	ID:stampy	VN:1.0.23_(r2059)	CL:--substitutionrate=0.01 -g /tmp/usecase3938954872417714817dir/references/R00000042/R00000042 -h /tmp/usecase3938954872417714817dir/references/R00000042/R00000042 -M bam -o /tmp/usecase3938954872417714817dir/0564a575-a6f5-40bc-8898-b0b5e944c4d9.sam --logfile=/tmp/usecase3938954872417714817dir/0564a575-a6f5-40bc-8898-b0b5e944c4d9.sam.log --readgroup=ID:824f45e8-37f3-4cb9-8a05-63f0b7c9b959 --outputformat=sam -v 3
@CO	PN:stampy	ID:stampy	VN:1.0.23_(r2059)	CL:--substitutionrate=0.01 -g /tmp/usecase3938954872417714817dir/references/R00000042/R00000042 -h /tmp/usecase3938954872417714817dir/references/R00000042/R00000042 -M bam -o /tmp/usecase3938954872417714817dir/570a7dde-6c04-419c-a898-87f872dd4eda.sam --logfile=/tmp/usecase3938954872417714817dir/570a7dde-6c04-419c-a898-87f872dd4eda.sam.log --readgroup=ID:7f568ff7-e0f6-4a55-ad17-6fe778ed8f83 --outputformat=sam -v 3
@CO	PN:stampy	ID:stampy	VN:1.0.23_(r2059)	CL:--substitutionrate=0.01 -t 8 -g /tmp/R00000042 -h /tmp/R00000042 --readgroup=ID:WTCHG_246141_245101,SM:7c2f06_45,PL:ILLUMINA,PU:160129_D00417_0381_AHJ2VGBCXX_1,LB:VAU2662A45,DT:2016-01-31,CN:WTCHG --comment=@MISC/WTCHG_246141_245101.comments.txt -M FASTQ/WTCHG_246141_245101_1.fastq.gz,FASTQ/WTCHG_246141_245101_2.fastq.gz
@CO	ID:stampy	TM:Sun, 31 Jan 2016 13:04:28 GMT	WD:/data1/GA-DATA/160129_D00417_0381_AHJ2VGBCXX/Data/Intensities/BaseCallsHN:comp03.mgmt.cluster2	UN:johnb
@CO	PN:stampy	ID:stampy	VN:1.0.23_(r2059)	CL:--substitutionrate=0.01 -t 8 -g /tmp/R00000042 -h /tmp/R00000042 --readgroup=ID:WTCHG_246142_245101,SM:7c2f06_45,PL:ILLUMINA,PU:160129_D00417_0381_AHJ2VGBCXX_2,LB:VAU2662A45,DT:2016-01-31,CN:WTCHG --comment=@MISC/WTCHG_246142_245101.comments.txt -M FASTQ/WTCHG_246142_245101_1.fastq.gz,FASTQ/WTCHG_246142_245101_2.fastq.gz
@CO	ID:stampy	TM:Mon, 07 Mar 2016 21:33:54 GMT	WD:/tmp/usecase3938954872417714817dir	HN:gel-pipeline3	UN:compass
@CO	PN:stampy	ID:stampy	VN:1.0.23_(r2059)	CL:--substitutionrate=0.01 -g /tmp/usecase3938954872417714817dir/references/R00000042/R00000042 -h /tmp/usecase3938954872417714817dir/references/R00000042/R00000042 -M bam -o /tmp/usecase3938954872417714817dir/a00a7733-2cfc-46b3-a685-3657fdee6848.sam --logfile=/tmp/usecase3938954872417714817dir/a00a7733-2cfc-46b3-a685-3657fdee6848.sam.log --readgroup=ID:824f45e8-37f3-4cb9-8a05-63f0b7c9b959 --outputformat=sam -v 3
@CO	ID:stampy	TM:Sun, 31 Jan 2016 15:02:32 GMT	WD:/data1/GA-DATA/160129_D00417_0381_AHJ2VGBCXX/Data/Intensities/BaseCallsHN:comp01.mgmt.cluster2	UN:johnb
@CO	CMD:/home/compass/PIPELINE/mmmPipeline/compass/g4_stampy.py -b bam -r R00000042 -o output -ss seqstats -fs flagstats -g e865d957-12e5-479a-9a08-131dfa0e9a5e""")

reads_string = """HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	99	R00000042	1619836	99	151M	=	1620303	618	CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC	DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII	PQ:i:205	SM:i:96	UQ:i:78	MQ:i:96	XQ:i:270	NM:i:2	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	99	R00000042	1619836	99	151M	=	1620303	618	CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC	DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII	PQ:i:205	SM:i:96	UQ:i:78	MQ:i:96	XQ:i:270	NM:i:2	RG:Z:test-test
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	1123	R00000042	1619836	99	151M	=	1620303	618	CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC	DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII	PQ:i:205	SM:i:96	UQ:i:78	MQ:i:96	XQ:i:270	NM:i:2	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	1123	R00000042	1619836	99	151M	=	1620303	618	CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC	DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII	PQ:i:205	SM:i:96	UQ:i:78	MQ:i:96	XQ:i:270	NM:i:2	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	147	R00000042	1620303	99	151M	=	1619836	-618	TATGAACGTCAGCTTTTGTGGTGATAAAGCTGGTGCGACGCTGTTTCCAGCACATATCTACACCGATAACAACGGTGAATTAATGACGCTGATGCAACGGGAAATGGCAGACGACAACCGCCATTTGCGCAGCATGGAAGCCTTTATCAAT	.HEF@7.AEF@@.@HG?EHHGCEEIHHEIIHIIHHIHHDHC@CFHIH@F70HEHEHCIIHHGHHHEIIIIIHIHGIIIIIIIIIIIHHHIIIIIFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD	PQ:i:205	SM:i:96	UQ:i:270	MQ:i:96	XQ:i:78	NM:i:7	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	147	R00000042	1620303	99	151M	=	1619836	-618	TATGAACGTCAGCTTTTGTGGTGATAAAGCTGGTGCGACGCTGTTTCCAGCACATATCTACACCGATAACAACGGTGAATTAATGACGCTGATGCAACGGGAAATGGCAGACGACAACCGCCATTTGCGCAGCATGGAAGCCTTTATCAAT	.HEF@7.AEF@@.@HG?EHHGCEEIHHEIIHIIHHIHHDHC@CFHIH@F70HEHEHCIIHHGHHHEIIIIIHIHGIIIIIIIIIIIHHHIIIIIFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD	PQ:i:205	SM:i:96	UQ:i:270	MQ:i:96	XQ:i:78	NM:i:7	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	1171	R00000042	1620303	99	151M	=	1619836	-618	TATGAACGTCAGCTTTTGTGGTGATAAAGCTGGTGCGACGCTGTTTCCAGCACATATCTACACCGATAACAACGGTGAATTAATGACGCTGATGCAACGGGAAATGGCAGACGACAACCGCCATTTGCGCAGCATGGAAGCCTTTATCAAT	.HEF@7.AEF@@.@HG?EHHGCEEIHHEIIHIIHHIHHDHC@CFHIH@F70HEHEHCIIHHGHHHEIIIIIHIHGIIIIIIIIIIIHHHIIIIIFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD	PQ:i:205	SM:i:96	UQ:i:270	MQ:i:96	XQ:i:78	NM:i:7	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	1171	R00000042	1620303	99	151M	=	1619836	-618	TATGAACGTCAGCTTTTGTGGTGATAAAGCTGGTGCGACGCTGTTTCCAGCACATATCTACACCGATAACAACGGTGAATTAATGACGCTGATGCAACGGGAAATGGCAGACGACAACCGCCATTTGCGCAGCATGGAAGCCTTTATCAAT	.HEF@7.AEF@@.@HG?EHHGCEEIHHEIIHIIHHIHHDHC@CFHIH@F70HEHEHCIIHHGHHHEIIIIIHIHGIIIIIIIIIIIHHHIIIIIFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD	PQ:i:205	SM:i:96	UQ:i:270	MQ:i:96	XQ:i:78	NM:i:7	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
HISEQ2500-09:381:HJ2VGBCXX:2:2201:15073:80781	1121	R00000042	1	70	1M7I143M	=	47901254790275	ATTTTTCAGCTTTTCATTCTGACTGCAATGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTCTCTGACAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATA	DDDDDHEHIHIIIIGHIIIIHHIIIIGIIHHFHHIHIFECFHHIGGHIIIFHEHIIIIIIIIHIHEHHIIIIIHHIIHIGIGH?HEHHHIIFGHHFHHHIHEFHIIIIIIHIIIIGHHHHHEHFHHHHHHHHHHEHHGHIFHHIGCHGHHH	PQ:i:375	SM:i:70	UQ:i:217	MQ:i:96XQ:i:186	NM:i:10	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959"""


def sam_record_from_sring(string):
    return pysam.AlignedSegment().fromstring(string, HEADER)


class TestGetReadIds(unittest.TestCase):

    def setUp(self):
        self.reads = []
        for line in reads_string.split('\n'):
            self.reads.append(sam_record_from_sring(line.strip()))

    def test_TwoReadIds_OnlyOneUnique(self):
        result = bam2fastq.get_read_ids(self.reads[:2])
        expected = {
            "HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635"
        }
        self.assertSetEqual(result, expected)

    def test_EightReadIds_OnlyTwoUnique(self):
        result = bam2fastq.get_read_ids(self.reads)
        expected = {
            "HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635",
            "HISEQ2500-09:381:HJ2VGBCXX:2:2201:15073:80781"
        }
        self.assertSetEqual(result, expected)


class TestGetUniqueReadPairs(unittest.TestCase):

    def setUp(self):
        self.reads = []
        for line in reads_string.split('\n'):
            self.reads.append(sam_record_from_sring(line.strip()))

        for i in range(len(self.reads)):
            self.reads[i].set_tag('RG', '824f45e8-37f3-4cb9-8a05-63f0b7c9b959')

    def test_OneR1_ReturnOneReadAndEmptyList(self):
        r1, r2 = bam2fastq.get_unique_reads_pairs(self.reads[:1])
        expected_r1 = [self.reads[0]]
        expected_r2 = []

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)

    def test_OneR2_ReturnOneReadAndEmptyList(self):
        r1, r2 = bam2fastq.get_unique_reads_pairs(self.reads[4:5])
        expected_r1 = []
        expected_r2 = [self.reads[4]]

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)

    def test_EmptyReads_ReturnEmptyLists(self):
        r1, r2 = bam2fastq.get_unique_reads_pairs([])
        expected_r1 = []
        expected_r2 = []

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)

    def test_OneR1OneR2_ReturnOneR1AndOneR2(self):
        r1, r2 = bam2fastq.get_unique_reads_pairs(
            [self.reads[0], self.reads[4]])
        expected_r1 = [self.reads[0]]
        expected_r2 = [self.reads[4]]

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)

    def test_MultipleDuplicates_ReturnOneR1AndOneR2(self):
        r1, r2 = bam2fastq.get_unique_reads_pairs(self.reads[:9])
        expected_r1 = [self.reads[0]]
        expected_r2 = [self.reads[4]]

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)

    def test_MultipleDuplicatesNoPrimary_ReturnOneR1AndOneR2(self):
        r1, r2 = bam2fastq.get_unique_reads_pairs(
            [self.reads[2], self.reads[3], self.reads[6]])
        expected_r1 = [self.reads[2]]
        expected_r2 = [self.reads[6]]

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)

    def test_MultipleDuplicatesTwoReadGroups_ReturnThreeReads(self):
        self.reads[1].set_tag('RG', 'test')
        r1, r2 = bam2fastq.get_unique_reads_pairs(self.reads)
        expected_r1 = self.reads[:2]
        expected_r2 = [self.reads[4]]

        self.assertListEqual(r1, expected_r1)
        self.assertListEqual(r2, expected_r2)


class TestToFastqString(unittest.TestCase):

    def test_EmptyRecord_EmptyString(self):
        result = bam2fastq.to_fastq_string(pysam.AlignedSegment())
        expected = ""

        self.assertEqual(result, expected)

    def test_FullRecord_CorrectFastqString(self):
        read = sam_record_from_sring(reads_string.split('\n')[0].strip())

        result = bam2fastq.to_fastq_string(read)
        expected = """@HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635/1 RG:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC
+
DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII"""

        self.assertEqual(result, expected)

    def test_RecordWithoutRG_CorrectFastqString(self):
        read = sam_record_from_sring("HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	99	R00000042	1619836	99	151M	=	1620303	618	CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC	DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII	PQ:i:205	SM:i:96	UQ:i:78	MQ:i:96	XQ:i:270	NM:i:2")

        result = bam2fastq.to_fastq_string(read)
        expected = """@HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635/1 
CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC
+
DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII"""

        self.assertEqual(result, expected)

    def test_RecordWithoutAssignedReadPair_CorrectFastqString(self):
        read = sam_record_from_sring("HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635	35	R00000042	1619836	99	151M	=	1620303	618	CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC	DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII	PQ:i:205	SM:i:96	UQ:i:78	MQ:i:96	XQ:i:270	NM:i:2	RG:Z:824f45e8-37f3-4cb9-8a05-63f0b7c9b959")

        result = bam2fastq.to_fastq_string(read)
        expected = """@HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635 RG:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC
+
DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII"""

        self.assertEqual(result, expected)


class TestWriteReadsToFastq(unittest.TestCase):

    def setUp(self):
        self.reads = []
        for line in reads_string.split('\n'):
            self.reads.append(sam_record_from_sring(line.strip()))

    def test_TwoUniquwRecords_TwoFastqEntries(self):
        outfile = StringIO()
        bam2fastq.write_reads_to_fastq([self.reads[0], self.reads[4]], outfile)

        outfile.seek(0)
        result = outfile.read()
        expected = """@HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635/1 RG:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
CCAGAACAGGCGCGGGAAATGTGCGATACCGCGCGCAAACTGGGCAAGGTGCTGGCCTACGACTTTCACCATCGTTTTGCGCTCGATACGCAACAGCTGCGTGAACAGGTGACCAACGGCGTTTTGGGAGAGATTTACGTTACCACCGCCC
+
DDDDDIIIIIIIIIIIIIIIHIHIIIIIIIIGIIIIIIIIGIIIHHIIIIGIIIIIIIIIIHIIIIIIIIIIIIIIIIIHIICGHIIIIHGIIIIHIIIGIGHIIIIIIIIIHHIIIIGIICHHHHIIHEHIIIIIIIIHHIIIIIHHIII
@HISEQ2500-09:381:HJ2VGBCXX:2:1101:10005:7635/2 RG:824f45e8-37f3-4cb9-8a05-63f0b7c9b959
TATGAACGTCAGCTTTTGTGGTGATAAAGCTGGTGCGACGCTGTTTCCAGCACATATCTACACCGATAACAACGGTGAATTAATGACGCTGATGCAACGGGAAATGGCAGACGACAACCGCCATTTGCGCAGCATGGAAGCCTTTATCAAT
+
.HEF@7.AEF@@.@HG?EHHGCEEIHHEIIHIIHHIHHDHC@CFHIH@F70HEHEHCIIHHGHHHEIIIIIHIHGIIIIIIIIIIIHHHIIIIIFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD
"""
        self.assertEqual(result, expected)

    def test_NoRecords_NothingWritten(self):
        outfile = StringIO()
        bam2fastq.write_reads_to_fastq([], outfile)

        outfile.seek(0)
        result = outfile.read()
        expected = ''
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
