from pathlib import Path

from .genomic_sequences_utils import fetch_seq

__all__ = ['Cfetch']

class Cfetch:
    """Chromosome data fetching."""

    def __init__(self, species: str = 'homo_sapiens') -> None:
        self.species = species
        self._rest_assembly: str = 'GRCh38'  # this should not be set by the user

    def fetch(self, chrm_info: Path | str, start: int, end: int, rev: bool = False) -> str:
        """Fetch sequence either from a local Fasta file or from Ensembl REST API.

        Args:
            chrm_info (Path | str): a chromosome fasta file or the chromosome number (e.g., '1' or 'Y').
            start (int): 1-based start coordinate
            end (int): 1-based end coordinate
            rev (bool, optional): reverse-complement the sequence. Defaults to False.

        Returns:
            str: the fetched sequence
        """
        return fetch_seq(chrm_info, start, end, rev=rev, species=self.species, rest_assembly=self._rest_assembly)
