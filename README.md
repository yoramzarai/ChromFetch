# Vertebrates Chromosome Sequence Retrieval

A Python tool to fetch an arbitrary sequence from a vertebrate chromosome.

# Python

Requires python>=3.12

# Details
The tool can fetch the sequence from either:
- a provided chromosome Fasta file, or
- remotely from Ensembl using the Ensembl REST API

In case of a provided chromosome Fasta file, the code loads only the queried sequence and not the complete Fasta file.

# Usage

```python
from pathlib import Path
from utils import Cfetch

# homo sapiens
cfetch = Cfetch(species='homo_sapiens')

# fetch from a provided chromosome Fasta file, by providing a chromosome Fasta file of type Path
fasta_file: Path = Path("./Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa")
seq = cfetch.fetch(fasta_file, 122_989_200, 122_989_200 + 29)
# or fetch from Ensembl, by providing a chromosome name of type str (e.g., '1', or 'Y')
seq = cfetch.fetch('1', 122_989_200, 122_989_200 + 29)

# mus musculus
cfetch = Cfetch(species='mus_musculus')
seq = cfetch.fetch('8', 11_100_100, 11_100_150, rev=True)
```

See also `ChromFetch.py`.