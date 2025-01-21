import os
from pathlib import Path

import pytest

from utils import Cfetch

# chromosomes main path on my local machine
local_human_chrm_path: Path = Path('/Users/yoramzarai/work/Databases/Reference_Genome/GRCh38.p14/Chromosome')
local_mouse_chrm_path: Path = Path('/Users/yoramzarai/work/Databases/Reference_Genome/Mus_musculus/GRCm39/Chromosome')

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

@pytest.mark.skipif(os.uname().nodename != 'Yorams-MacBook-Pro-M1-MAX.local', reason="Local Fasta files only on my MAC.")
@pytest.mark.parametrize(
    "fasta_file, start, end, species, rev, expected",
    [
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa', 123_123_123, 123_123_233, 'homo_sapiens', False, "TGGTAGAAAAGGAAATCTCTTCGTATAAAAACTAGACAGAATCACTCTCAGAAACTGCTCTGCGATGTGTGCGTTCAACTCTCAGAGTTTAACTTTTCTTTTCATTCAGCA"),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa', 23_4567_890, 23_4567_990, 'homo_sapiens', True, 'ACTTTTCCCTGTCCCAATAAGGGCATCTTATGTCTGGGTGTCATTTTCATCAGAATGCAAAGGTAAACAAGTGAGAAGCGTCTTTGTGTGGCTGGAGCTTA'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.3.fa', 45_100_100, 45_100_210, 'homo_sapiens', True, 'TTTCTTCCCTTCGGTAGAGAATGACTTCTTCTCTGGGGGCATCTGACCAATCTAAGAGAAAAGACCCAAAGATATTGATGTTATGGGCTCCCCAGGAAGCAATCTGGCCAT'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.4.fa', 876543, 876643, 'homo_sapiens', False, 'CAGTGGTGGCGTTCCGAGGCACAAACCCCGTGTGGAACTGAATCTGGAACATCTTCATGGATGCCATCTGCAAAGAGAGCAAACACGACACCCCACGTGGA'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.5.fa', 1_900_876, 1_900_976, 'homo_sapiens', True, 'ACTCAAGGAATGAGCAACAGAAGGAGAAGACAAGCCAGCAAAGGGCTCCTGCATCCCTGGGTTTTAAGTCTGTTTCCCCGGCTCTAGGTCTCGGAGAGCTC'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.6.fa', 99_987_100, 99_987_200, 'homo_sapiens', False, 'ATTTATTTATTTATTTTTTGAGACAGAGTCTTGCTCTGCTGCCCAGGGTGGAGTGCAGCGGTGTGATCTCGGCTCACTGCAACCTCTGCCTCCTGGGTTCA'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa', 100999, 101100, 'homo_sapiens', True, 'CAGGGTGTTGGATTTTCTCAAGTCCTTTTTCTCTATCAATTGAGACAGTCGTGTAGTTTTTGTCCTTCATTATGACAGTGTGGCAGATTATCTTGATTGATG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.8.fa', 455_911, 455_999, 'homo_sapiens', False, 'ATAGCCATTAACAGGGATGTAAAACTGGACTCTGCCAGAAATCCCAAATGCTGAGTTGTCACACCAGACTCTATAACTCTACTAATCCG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.9.fa', 29_876, 29_991, 'homo_sapiens', False, 'CACAACGGTGTGAATCTACTTAATCCCACTGAACTGTATGCTGAAAAATGGTTTAGACGGTGAATTTTAGGTTATGTATGTTTTACCACAATTTTTAAAAAGCTAGTGAAAAGCTG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.10.fa', 100998, 101091, 'homo_sapiens', True, 'CGTATTTTTTTGTAGATATGAGGTTTCACCATGTTTCCCAGGCTGGTCTCAAACTCCCCTGCTCAAGCAATATTAACACCTTGGCCTCCCAAAG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.11.fa', 60_000_000, 60_000_100, 'homo_sapiens', False, 'TTCTGAGGCAGTCTCCCTGGCTGGAAGCAGTCCTTACTGTTAAGGACTCACTCAATTAAATTGGGTTTCTCAGACAATCTGGGATAATCTCTCTATTTTAA'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.12.fa', 60_000_001, 60_000_100, 'homo_sapiens', False, 'AGGAATCACCAATGTTCCAAACAAAATGTCATAAAAAGACATTTTATTCAAATGTAATGAATTTTTGTTTTTTGGAATTGAGTTGGCAGTTTCATGAACC'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.13.fa', 59_999_999, 60_000_100, 'homo_sapiens', False, 'TAACAGAACACTTTTGCTTTTTCTTCTTGTTCTTTATATATGACCAGCTTCTGTATTTTCATATGTTACTATGAGATACGGGTTAAGACCATGTGTATGCCA'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.14.fa', 59_999_989, 60_000_020, 'homo_sapiens', False, 'ATTAATATGTTCATTGCAGCACTATTCACAAT'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.15.fa', 90_000_000, 90_000_060, 'homo_sapiens', False, 'CCCTTGCATTGTTCCCTGAGTTTGCTCTCGGGGCCCACCGAGACCTCACGTTTCTCCCCCT'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.16.fa', 90_000_001, 90_000_060, 'homo_sapiens', True, 'CTCCTGCGGCCCGTGCCGGGGGATGCCCACCCTTACCTGGGCGTAGGCACTCTGGGTGAC'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa', 60_000_000, 60_000_170, 'homo_sapiens', True, 'CATTTTAGGCCAGTTTGAAATAGAAAAAAGTTGTGAGTTACTTTTTCAGATGGCAGAGTCTGTTGATTATGATTATGATGTTCAAGATACCATCCTCGATGGCTTGTTAATATTGCCGTGTTATGGATCAATGACAACTGGTAATTTCTCATTAGAATAGAAAATTTGATT'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.18.fa', 60_000_000-10, 60_000_120, 'homo_sapiens', True, 'AAACCAGGAAGAAGCTGAATTCCTGAATAGATCAATAACTAGTTGTGAAACTGAGGCAGTAATTAATAGCCTACCAACCAAAAAATAAGCCCAGGACCAAACAGATTCACAGCCAAATTCTACCAGAGGTA'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.19.fa', 123461, 123499, 'homo_sapiens', False, 'AGCCAAAAAACAGATCAAATCAGTAAACCAAAAATCTTG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa', 60002, 60061, 'homo_sapiens', False, 'GTTCAGTCGGGCAGGGAGTGGGAATAGACAAGACCACAAGCAGCTTGGTGCCTCTGAAAG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.21.fa', 37365981, 37365999, 'homo_sapiens', False, 'TCAGGGCCGCCGGGCCCGG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.22.fa', 37365981, 37365999, 'homo_sapiens', True, 'CCCCAATTCTGAATTAATG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.X.fa', 99_879, 99_987, 'homo_sapiens', True, 'CACATAGAGAGCAGACTGCGCAAACTTTAGAGTCTGCATTGGGCCTAGGTCTCACTGAGGACAGATAGAGAGCACACTGTGTAACCTTTAGAGTCTGCATAGGGCCTCG'),
        (local_human_chrm_path / 'Homo_sapiens.GRCh38.dna_sm.chromosome.Y.fa', 123409, 123511, 'homo_sapiens', False, 'CCGCAATTAGACCTAGGCCCAATGCAGACTCTAAAGGTTGCAGAGTCTGCTGTCTATCTGACCTCAAGGAGACCTAGGCCCAATGCAGACTCTAAAGGTTGCA'),

        # mouse
        (local_mouse_chrm_path / 'Mus_musculus.GRCm39.dna_sm.chromosome.2.fa.gz', 60_000_000, 60_000_123, 'homo_sapiens', False, 'GAATGTCTCCTGCCTTGGCCTCCCTTGTGCTGAGAAATCTATGGTTAGTATTTGTCTTAGTCACTGTTTTACTGTGTGAAGAGACACCATGGCCAAGGTATGTCTTATAAAAGGAAGCATTTAA'),
        (local_mouse_chrm_path / 'Mus_musculus.GRCm39.dna_sm.chromosome.4.fa.gz', 136_366_473, 136_366_573, 'homo_sapiens', True, 'TAAAAAGCAATGATGAACTATGGACACACACAAAGACTCGGACCTCAAAACCAGCACGCTGAGCCCAAGGAGATGTACCTAAAACCTCTGTGTGCTCCCAT'),
    ],
)
def test_with_local_fasta(fasta_file: Path, start: int, end: int, species: str, rev: bool, expected: str) -> None:
    seq = Cfetch(species=species).fetch(fasta_file, start, end, rev=rev)
    assert seq.upper() == expected
