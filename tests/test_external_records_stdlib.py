"""Bounded, explicitly CC0 external source parsing; no competition inputs needed."""
import unittest

from casmi.validate.massbank_records import parse_record, RecordRejected


RECORD = """ACCESSION: EXTERNAL-SYNTHETIC-1
AUTHORS: Synthetic test fixture
LICENSE: CC0
CH$SMILES: CCO
AC$INSTRUMENT_TYPE: LC-ESI-QTOF
AC$MASS_SPECTROMETRY: MS_TYPE MS2
AC$MASS_SPECTROMETRY: IONIZATION ESI
AC$MASS_SPECTROMETRY: ION_MODE POSITIVE
MS$FOCUSED_ION: PRECURSOR_M/Z 47.0491
MS$FOCUSED_ION: PRECURSOR_TYPE [M+H]+
PK$NUM_PEAK: 5
PK$PEAK: m/z int. rel.int.
  10 1 100
  11 2 200
  12 5 500
  13 10 999
  14 1 100
//
"""


class ExternalRecordTests(unittest.TestCase):
    def test_cc0_attribution_and_intensity_normalization(self):
        row = parse_record(RECORD)
        self.assertEqual(row["license"], "CC0")
        self.assertEqual(row["authors"], "Synthetic test fixture")
        self.assertEqual(row["ms2_normalized_intensities"], [10, 20, 50, 100, 10])
        self.assertEqual(len(row["record_sha256"]), 64)

    def test_ambiguous_unlicensed_or_invalid_records_fail_closed(self):
        mutations = [RECORD.replace("LICENSE: CC0", "LICENSE: CC BY"),
                     RECORD.replace("LICENSE: CC0", "LICENSE: CC0\nLICENSE: CC0"),
                     RECORD.replace("MS_TYPE MS2", "MS_TYPE MS3"),
                     RECORD.replace("IONIZATION ESI", "IONIZATION EI"),
                     RECORD.replace("PK$NUM_PEAK: 5", "PK$NUM_PEAK: 6"),
                     RECORD.replace("  12 5 500", "  12 nan 500"),
                     RECORD.replace("  12 5 500", "  11 5 500"),
                     RECORD.replace("47.0491", "nan"), RECORD.replace("//", "")]
        for index, text in enumerate(mutations):
            with self.subTest(case=index), self.assertRaises(RecordRejected):
                parse_record(text)


if __name__ == "__main__":
    unittest.main()
