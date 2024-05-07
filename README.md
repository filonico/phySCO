# Welcome to phySCO!

phySCO is a python script that infers maximum likelihood <ins>phy</ins>logenomic tree using BUSCO <ins>s</ins>ingle-<ins>c</ins>opy <ins>o</ins>rthologous genes. phySCO retrieves these genes from already-available BUSCO results (that is, phySCO *does not* run BUSCO itself).

Many thanks to [NiccoloRighetti](https://github.com/NiccoloRighetti), whose work has pushed me to start writing this code.

## General description
Given a directory containing BUSCO results for a group of species, phySCO computes the maximum likelihood phylogenetic tree of all the species.

To run phySCO, just create a dedicated directory and copy into it BUSCO results for the species you want to analyse.
> **Mind to keep the default structure of BUSCO results**.

phySCO is currently able to process **just amino acid sequences** of complete single-copy BUSCO genes. In fact, amino acids sequences are always returned by any BUSCO analysis run mode (either genome, transcriptome or protein). Maybe, in the future, phySCO will be able to process also nucleotide sequences...

## Required softwares and dependencies
Here is the list of software that phySCO requires:
* <code>python</code> (v3.11)
* <code>mafft</code> (v7.505)
* <code>trimal</code> (v1.4.rev15)
* <code>IQTREE2</code> (2.1.4-beta COVID-edition)

You can install all of them directly from the [<code>phySCO_env.yml</code>](./phySCO_env.yml) YAML file, through the command:
```
conda env create -f phySCO_env.yml
conda activate phySCO_env
```

## Example dataset
[<code>example_dataset/</code>](./example_dataset/) is a test dataset that can be used directly as an input to phySCO. [<code>example_dataset_key.md</code>](./example_dataset_key.md) contains metadata of the example dataset.

It has been generated by simply running BUSCO on a random set of mammal NCBI reference genomes. You can find the keys to the species identifier in [<code>example_dataset_key.md</code>](./example_dataset_key.md).

Before running phySCO on the example dataset, run the following commands:
```
mkdir example_dataset_extracted
for i in example_dataset/*tar.gz; do tar -xvzf $i -C example_dataset_extracted; done
```
then run
```
python3 phySCO.py -i example_dataset_extracted/
```
