source $CONDA_PREFIX/etc/profile.d/mamba.sh
source $CONDA_PREFIX/etc/profile.d/conda.sh

## Bnapus ZS11 setup

GFF=../validation/ZS11/zs11.genome.noScaf.gff
FASTA=../validation/ZS11/zs11.genome.fa
./prep --gff=$GFF --fasta=$FASTA --outdir=../validation/ZS11/data

## Bnapus ZS11 predictions

CUDA_VISIBLE_DEVICES=0 python nemo/predict.py ../validation/ZS11/data


## Athaliana setup

VAL=../validation/Athaliana
EQTL=${VAL}/data/eQTL
GATC=${VAL}/data/GATC
mkdir -p $EQTL
mkdir -p $GATC

GFF=../data/Athaliana/parent_data/TAIR10_GFF3_genes.gff


LIST=${VAL}/genes.eQTL.txt
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$EQTL"/"$STEM
  mkdir $OUTDIR
  echo $FASTA
  echo $OUTDIR
done | parallel -N 2 -j 10 ./prep --gff=$GFF --fasta={1} --outdir={2} --genelist=$LIST


# check if everythin was set up correctly
NUM_EXPECTED=$(wc -l $LIST | cut -f1 -d" ")
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$EQTL"/"$STEM
  if [[ $(wc -l $EQTL"/"$STEM/gene.id.key | cut -f1 -d" ") -ne $NUM_EXPECTED ]]; then
    rm -r $OUTDIR
    ./prep --gff=$GFF --fasta=$FASTA --outdir=$OUTDIR --genelist=$LIST
  fi
done

LIST=${VAL}/genes.GATC.txt
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$GATC"/"$STEM
  mkdir $OUTDIR
  echo $FASTA
  echo $OUTDIR
done | parallel -N 2 -j 10 ./prep --gff=$GFF --fasta={1} --outdir={2} --genelist=$LIST

# check if everythin was set up correctly
NUM_EXPECTED=$(wc -l $LIST | cut -f1 -d" ")
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$GATC"/"$STEM
  if [[ $(wc -l $GATC"/"$STEM/gene.id.key | cut -f1 -d" ") -ne $NUM_EXPECTED ]]; then
    rm -r $OUTDIR
    ./prep --gff=$GFF --fasta=$FASTA --outdir=$OUTDIR --genelist=$LIST
  fi
done


#### Athaliana predictions

export count=0; for file in ../validation/Athaliana/*.fa; do
  if [[ count -le 316 ]]; then
    echo $file >> ../validation/Athaliana/genomes_set1.txt
    let count++
  else
    echo $file >> ../validation/Athaliana/genomes_set2.txt
    let count++
  fi
done

modelfile=../model_weights/nemo/Athaliana/masked_graphpart/nemo_full.h5

mamba activate nemo
for fasta in $(cat ../validation/Athaliana/genomes_set1.txt)
do
    name=$(basename $fasta .chr.fa)
    dir=../validation/Athaliana/data/eQTL/$name
    CUDA_VISIBLE_DEVICES=0 python nemo/predict.py $dir $modelfile
    dir=../validation/Athaliana/data/GATC/$name
    CUDA_VISIBLE_DEVICES=0 python nemo/predict.py $dir $modelfile
done & pid1=$!

for fasta in $(cat ../validation/Athaliana/genomes_set2.txt)
do
    name=$(basename $fasta .chr.fa)
    dir=../validation/Athaliana/data/eQTL/$name
    CUDA_VISIBLE_DEVICES=1 python nemo/predict.py $dir $modelfile
    dir=../validation/Athaliana/data/GATC/$name
    CUDA_VISIBLE_DEVICES=1 python nemo/predict.py $dir $modelfile
done & pid2=$!

wait $pid1 $pid2
mamba deactivate


## Athaliana aggregate
OUTFILE=../validation/Athaliana/validation.eQTL.tsv
echo -e "accession_id\tgene_id\tpredicted_expression" > $OUTFILE
for DIR in ../validation/Athaliana/data/eQTL/*
do
  genotype=$(echo $DIR | rev | cut -f1 -d"/" | rev)
  FILE=$DIR"/predictions.txt"
  awk -v genotype="$genotype" 'BEGIN{OFS="\t"};NR>1{print genotype, $1, $2}' $FILE >> $OUTFILE
done

OUTFILE=../validation/Athaliana/validation.GATC.tsv
echo -e "accession_id\tgene_id\tpredicted_expression" > $OUTFILE
for DIR in ../validation/Athaliana/data/GATC/*
do
  genotype=$(echo $DIR | rev | cut -f1 -d"/" | rev)
  FILE=$DIR"/predictions.txt"
  awk -v genotype="$genotype" 'BEGIN{OFS="\t"};NR>1{print genotype, $1, $2}' $FILE >> $OUTFILE
done
