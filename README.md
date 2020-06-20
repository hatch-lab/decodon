# Decodon

## Installation
```
git clone https://github.com/hatch-lab/decodon.git
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Usage
`python decodon.py "DNA Sequence"`

This will print to terminal DNA sequences that are as divergent as possible, but which encode the same amino acid sequence.

#### Options
- `--N` The number of divergent sequences to print (defaults to 1)

#### Examples
```
python decodon.py GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG --N=3
```
Prints
```
Sequence (Score)
GGCAAACCCATTCCCAATCCCTTATTAGGCTTAGACAGCACT (19)
GGAAAGCCAATACCAAACCCATTGTTGGGATTGGATAGTACC (15)
GGGAAACCGATCCCGAATCCGCTTCTTGGGCTTGACTCCACA (13)
```

## Output
Prints DNA sequences to the terminal.