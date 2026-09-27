The *napus* expression model (*n*emo) is a deep convolutional neural network that predicts median gene expression across tissues under control conditions in rapeseed (*Brassica napus*).

![nemo logo](logo/nemo.jpg?raw=true "The napus expression model")

# Installation

```
wget https://github.com/KRockenbach/nemo/nemo.tar.gz
tar -xzvf nemo.tar.gz
cd nemo
```

Environments are managed using Mamba (https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).
For GPU acceleration, nvidia drivers >=525.60.13 must be installed.
If no nvidia drivers are installed inference will run on CPU.
```
mamba env create -f nemo-predict.yml
```


# Make predictions on your own data

To make predictions, a genome assembly in fasta format and a structural annotation in gff3 format is needed.
The annotation should contain mRNA and CDS features. Each mRNA feature should have an ID tag in the attributes column.

```
mamba activate nemo-predict
python nemo.py -f/--fasta <ASSEMBLY> -g/gff <ANNOTATION> -p/--prefix <OUTPUT_PREFIX>
```

A tmp directory can be set using the optional `-t/--tmp` option, or by setting the TMPDIR environment variable.

The output is a tab separated file `<OUTPUT_PREFIX>.predicions.tsv` containing predicted TPM normalized expression values for each transcript.
Model outputs are rounded to 1 decimal and clipped at zero. 
The expression values represent median expression across Brassica napus tissues under standard greenhouse conditions.

Memory comsumption can be managed using the `-b/--batch_size` option. The default batch size is 120.


# How to cite us

Please cite as:

Kevin C Rockenbach, Silvia F Zanini, Alison C Tidy, Richard J Morris, Rachel Wells, Agnieszka A Golicz, A deep learning model captures position-specific effects of plant regulatory sequences and suggests genes under complex regulation, Plant Physiology, 2026;, kiag473, (https://doi.org/10.1093/plphys/kiag473)
