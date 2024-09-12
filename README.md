# ClustAlignPy

## Description
This program performs multiple sequence alignment (MSA) for proteins by combining dynamic programming (Needleman-Wunsch algorithm) with the UPGMA clustering method. The sequences are aligned progressively according to the order defined by a phylogenetic tree.

## Installation
To run the program, you need Python 3.x and the following dependencies:
- `numpy`
- `pandas`
- `biopython` (for handling FASTA files)

To install the dependencies, create a conda enviroment using the .yml file provided:

```bash
$ conda env create -f clustalignpy.yml
```

## Usage
Don't forget to activate the conda enviroment.
```bash
conda activate clust_align_py_env
```

You can run the program by passing a FASTA file containing the sequences you wish to align.

```bash
python multiple_alignment_NW.py sequences.fasta
```

## Steps

The program will:

1. Read the sequences from the FASTA file.
2. Perform pairwise alignments and construct a score matrix.
3. Convert the score matrix to a distance matrix and construct a phylogenetic tree using UPGMA.
4. Align the sequences progressively based on the UPGMA tree order.
5. Output the aligned sequences.

## Example Output
The output consists of the aligned sequences printed in a formatted manner in the console, showing gaps (-) and the number of amino acids aligned in each segment.

Sequence_1: M--KTADGV...    7
Sequence_2: MK-KTADGV...    8
Sequence_3: MK-LTAD--...    6


## Limitations
 - Currently, the program uses a fixed gap penalty. Incorporating affine gap penalties could provide more biologically realistic alignments.
 - It does not support parallelization, which might limit its use for very large datasets.

## Contact
For any issues or questions, feel free to reach out to [oumayma.meskini@etu.u-paris.fr](mailto:oumayma.meskini@etu.u-paris.fr).








