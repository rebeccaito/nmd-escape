import unittest
import os
from annotating_nmd import *


class TestNMD(unittest.TestCase):
    test_dir = os.path.dirname(os.path.abspath(__file__))
    test_bed_file = os.path.join(test_dir, 'data/test_cds.bed')

    cds_bed = pd.read_table(test_bed_file, names=['chrom', 'start', 'end', 'cds_id', 'score', 'strand'])
    cds_bed.sort_values('start', inplace=True)

    def test_preprocess_bed(self):
        out_columns = ['chrom', 'start', 'end', 'cds_id', 'score', 'strand', 'transcript_name', 'cds_size']

        # wrong number of columns, raise error
        with self.assertRaises(Exception) as context:
            preprocess_bed(TestNMD.cds_bed[['chrom', 'start']].copy())

        self.assertEqual(context.exception.args[0], "BED file is not in 6 column format")

        # different column names, rename
        with self.assertWarns(Warning) as context:
            renamed_cds_bed = TestNMD.cds_bed.rename(columns={'chrom': 'chr'}).copy()
            out_bed = preprocess_bed(renamed_cds_bed.copy())
            self.assertEqual(list(out_bed.columns), out_columns)
            self.assertTrue(issubclass(context.warnings[0].category, IncorrectColumnWarning))

        # runs correctly, adds size and transcript_name columns
        out_bed = preprocess_bed(TestNMD.cds_bed.copy())
        self.assertEqual(list(out_bed.columns), out_columns)

    def test_get_nmd_escape_size(self):
        test_df = preprocess_bed(TestNMD.cds_bed.copy())
        small_cds_df = test_df.copy()
        small_cds_df['cds_size'] = 40
        tests = [
            {
                'test_name': '1 CDS exon, + strand',
                'data': test_df[test_df.strand == '+'].iloc[[0]],
                'expected_size': 472
            }, {
                'test_name': '>1 CDS exon, + strand, penultimate exon >55',
                'data': test_df[test_df.strand == '+'],
                'expected_size': 558 + 55
            }, {
                'test_name': '>1 CDS exon, + strand, penultimate exon <55',
                'data': small_cds_df[small_cds_df.strand == '+'],
                'expected_size': 80
            },             {
                'test_name': '1 CDS exon, - strand',
                'data': test_df[test_df.strand == '-'].iloc[[0]],
                'expected_size': 2045
            }, {
                'test_name': '>1 CDS exon, - strand, penultimate exon >55',
                'data': test_df[test_df.strand == '-'],
                'expected_size': 2045 + 55
            }, {
                'test_name': '>1 CDS exon, - strand, penultimate exon <55',
                'data': small_cds_df[small_cds_df.strand == '-'],
                'expected_size': 80
            }, {
                'test_name': 'no CDS size provided (preprocess_bed is called)',
                'data': TestNMD.cds_bed[TestNMD.cds_bed.strand == '+'].copy(),
                'expected_size': 558 + 55
            }
        ]

        for test in tests:
            with self.subTest(msg=test['test_name']):
                result = get_nmd_escape_size(test['data'].copy())
                self.assertEqual(result, test['expected_size'])

    def test_make_cds_size_df(self):
        sizes_df = make_cds_size_df(TestNMD.cds_bed.copy())
        expected_cds_sizes = [1818, 2685]
        expected_nmd_sizes = [613, 2100]
        self.assertEqual(sizes_df.cds_size.to_list(), expected_cds_sizes)
        self.assertEqual(sizes_df.nmd_escape_size.to_list(), expected_nmd_sizes)

    def test_get_nmd_escape_boundaries(self):
        # test on 1 transcript
        nmd_df = get_nmd_escape_boundaries(TestNMD.cds_bed[TestNMD.cds_bed.strand == '+'].copy())
        self.assertEqual(len(nmd_df), 2)

        # test on all transcripts
        cds_df = preprocess_bed(TestNMD.cds_bed.copy())
        nmd_df = cds_df.groupby('transcript_name').apply(get_nmd_escape_boundaries)
        self.assertEqual(len(nmd_df), 4)

        # test convenience wrapper
        nmd_df = make_boundaries_df(TestNMD.cds_bed.copy())
        self.assertEqual(len(nmd_df), 4)

    def test_get_upstream_frameshift(self):

        test_variants = [
            {
                'test_name': 'upstream NMD escaping frameshift',
                'transcript_name': 'NM_003620.4',
                'HGVSp': 'NP_003611.1:p.Cys300LeufsTer129',
                'expected_result': True
            }, {
                'test_name': 'frameshift NMD inducing',
                'transcript_name': 'NM_003620.4',
                'HGVSp':'NP_003611.1:p.Cys81LeufsTer129',
                'expected_result': False
            }, {
                'test_name': 'frameshift NMD inducing',
                'transcript_name': 'NM_138576.4',
                'HGVSp':'NP_612808.1:p.Leu71TrpfsTer4',
                'expected_result': False
            }, {
                'test_name': 'upstream NMD escaping frameshift',
                'transcript_name': 'NM_138576.4',
                'HGVSp': 'NP_612808.1:p.Ala512ProfsTer21',
                'expected_result': True
            },{
                'test_name': 'not a frameshift (invalid)',
                'transcript_name': 'NM_138576.4',
                'HGVSp': 'NP_612808.1:p.Ala512Pro',
                'expected_result': False
            },{
                'test_name': 'improperly formatted frameshift',
                'transcript_name': 'NM_138576.4',
                'HGVSp': 'I AM A NAUGHTY VARIANT',
                'expected_result': False
            },{
                'test_name': 'frameshift not truncating',
                'transcript_name': 'NM_138576.4',
                'HGVSp': 'NP_612808.1:p.Ala512ProfsTer?',
                'expected_result': False
            }
        ]

        test_df = pd.DataFrame(test_variants)
        sizes_df = make_cds_size_df(TestNMD.cds_bed.copy())
        result = get_upstream_frameshift(test_df.copy(), sizes_df)
        self.assertTrue((result.expected_result == result.is_nmd_frameshift).all())
        self.assertEqual(len(result), 7)

        # throw warnings
        with self.assertWarns(Warning) as context:
            get_upstream_frameshift(test_df.iloc[-3:].copy(), sizes_df)
            self.assertTrue(any([issubclass(warning.category, HGVSpPatternWarning) for warning in context.warnings]))


if __name__ == '__main__':
    unittest.main()
