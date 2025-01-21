from pathlib import Path

from utils import Cfetch


def main() -> None:
    """Usage examples."""
    # homo sapiens
    cfetch = Cfetch(species='homo_sapiens')

    # fetch from a provided chromosome Fasta file, by providing a chromosome Fasta file of type Path
    fasta_file: Path = Path("/Users/yoramzarai/work/Databases/Reference_Genome/GRCh38.p14/Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa")
    seq = cfetch.fetch(fasta_file, 122_989_200, 122_989_200+30-1)
    print(f"{seq=}")

    # fetch from Ensembl, by providing a chromosome name of type str (e.g., '1', or 'Y')
    seq = cfetch.fetch('1', 122_989_200, 122_989_200+30-1)
    print(f"{seq=}")

    # Danio rerio
    cfetch = Cfetch(species='Danio_rerio')
    seq = cfetch.fetch('5', 52_120_100, 52_120_120, rev=True)
    print(f"{seq=}")

    # mus musculus
    cfetch = Cfetch(species='mus_musculus')
    seq = cfetch.fetch('8', 11_100_100, 11_100_150, rev=False)
    print(f"{seq=}")


if __name__ == "__main__":
    main()
