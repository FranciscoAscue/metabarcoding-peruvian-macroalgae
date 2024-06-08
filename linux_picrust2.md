#### Install Anaconda
#### Create Picrust2 environment

```bash
wget https://raw.githubusercontent.com/picrust/picrust2/master/picrust2-env.yaml
conda env create -f picrust2-env.yaml
conda activate picrust2
sudo apt install samtools
### cd ./results/picrust2/
# Picrust2 files
### metadata.tsv
### asv.table
### asv-seq.fna
### otu_norarefy.biom
```

#### Filter fasta for contaminants

```bash
sed -i 's/"//g' asv.table
samtools faidx asv-seq.fna
samtools faidx asv-seq.fna $(cat asv.table)
samtools faidx asv-seq.fna $(cat asv.table) > asv-seq_filter.fna
```

#### Add metadata to BIOM format and run Picrust2

```bash
head -n 2 metadata.tsv | sed "s/\\t/,/g" > names.txt
biom add-metadata -i otu_norarefy.biom -o otu_norarefy_metadata.biom -m metadata.tsv --sample-header $(cat names.txt)
picrust2_pipeline.py -s rep-seq_filter.fna -i table_with_sample_metadata.biom -o results -p 48
```