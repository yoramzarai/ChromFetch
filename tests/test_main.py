from pathlib import Path

import pytest

from utils import Cfetch


@pytest.mark.parametrize(
    "chrm, start, end, species, rev, expected",
    [
        ('1', 122_989_200, 122_989_200+29, 'homo_sapiens', False, "GATATTTTGACCACTTAGAGGCCTTCGTTG"),
        ('1', 112_989_200, 112_989_200+39, 'homo_sapiens', True, "TAAGTCTCATCCTAATTAAATCTTCTAGAATAAGCCTAAA"),
        ('5', 52_120_100, 52_120_120, 'Danio_rerio', True, "TTTTAGAGGTTGAACAGCCAC"),
    ],
)
def test_species_no_fasta(chrm: str, start: int, end: int, species: str, rev: bool, expected: str) -> None:
    seq = Cfetch(species=species).fetch(chrm, start, end, rev=rev)
    assert seq.upper() == expected


@pytest.mark.parametrize(
    "fasta_file, start, end, species, rev, expected",
    [
        (Path('./Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa'), 1_100_100, 1_100_150, 'homo_sapiens', False, "TTGAGACGGAGTTTTGGTCTTGTAGCCCAGACTGGAGTGCAATGGTGCAAT"),
        (Path('./Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa'), 1_200_100, 1_200_150, 'homo_sapiens', True, 'TTCACTTTACTTTGTTGGCTCACTCTTGAATTCTTTCCTTTATGAAGCCAA'),
    ],
)
def test_with_fasta(fasta_file: Path, start: int, end: int, species: str, rev: bool, expected: str) -> None:
    seq = Cfetch(species=species).fetch(fasta_file, start, end, rev=rev)
    assert seq.upper() == expected