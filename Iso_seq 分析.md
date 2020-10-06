

## Iso_seq 分析

### workdir

```shell
workdir=/share/home/stu_marui/feisha.test/iso_seq
```

### 1. pbindex

```shell 
pbindex $workdir/PacBio_Rawdata/m64012_181222_192540.subreads.bam
pbindex $workdir/PacBio_Rawdata/m64012_181221_231243.subreads.bam

```

### 2. ccs

```shell
#01_ccs
ccs \
--min-passes 3 \
--min-length 50 \
--max-length 15000 \
--min-rq0.99 \
$workdir/PacBio_Rawdata/m64012_181222_192540.subreads.bam \
$workdir/result/m64012_181222_192540.ccs.bam

ccs \
--min-passes 3 \
--min-length 50 \
--max-length 15000 \
--min-rq0.99 \ $workdir/PacBio_Rawdata/m64012_181221_231243.subreads.bam \ $workdir/result/m64012_181221_231243.ccs.bam


```

### 3. lima (CCS分类得到full-length reads)

```shell
#02_lima
lima \
--isoseq \
--dump-clips \
-j 8 \
$workdir/result/m64012_181222_192540.ccs.bam \
IsoSeqPrimers.fasta \
$workdir/result/m64012_181222_192540.fl.bam

lima \
--isoseq \
--dump-clips \
-j 8 \
$workdir/result/m64012_181221_231243.ccs.bam \
IsoSeqPrimers.fasta \
$workdir/result/m64012_181221_231243.fl.bam

#结果输出：
#m64012_181222_192540.fl.primer_5p--primer_3p.bam
#m64012_181221_231243.fl.primer_5p--primer_3p.bam

```

### 4. isoseq3_refine (去除polyA和嵌合序列)

```shell
#03_isoseq3_refine
isoseq3 refine \
--require-polya \
--min-polya-length 20 \
-j 4 \
$workdir/result/m64012_181222_192540.fl.primer_5p--primer_3p.bam \
IsoSeqPrimers.fasta \
$workdir/result/m64012_181222_192540.flnc.bam

isoseq3 refine \
--require-polya \
--min-polya-length 20 \
-j 4 \
$workdir/result/m64012_181221_231243.fl.primer_5p--primer_3p.bam \
IsoSeqPrimers.fasta \
$workdir/result/m64012_181221_231243.flnc.bam

```

### 5. 创建 flnc.xml

```shell
dataset create \
--type ConsensusReadSet \
$workdir/result/combined.flnc.xml \
$workdir/result/m64012_181222_192540.flnc.bam \
$workdir/result/m64012_181221_231243.flnc.bam
```

### 6. isoseq3 cluster (FLNC聚类)

```shell
isoseq3 cluster \
-j 8 \
--verbose \
--use-qvs \
$workdir/result/combined.flnc.xml \
$workdir/result/polished.bam 


#结果输出：
#polished.bam
#polished.lq.fasta.gz/polished.hq.fasta.gz
```



## 参考基因组比对

### mapping

##### 1) gmap 比对

```shell
reference= reference/path
workdir=

#01_index
gmap_build \
-D gmapdb \
-d hg38 \
-k 15 \
$reference/hg38/hg38.fa

#02_gmap
gmap \
-D gmapdb \
-d hg38 \
-f samse \
-n 0 \
-t 8 \
-z sense_force \
$workdir/data/flnc.fasta > $workdir/result/flnc.fasta.sam 2> flnc.fasta.sam.log

#03_sort
samtools view -bS $workdir/result/flnc.fasta.sam > $workdir/result/flnc.fasta.bam
samtools sort $workdir/result/flnc.fasta.bam > $workdir/result/flnc.fasta.sorted.bam
samtools index $workdir/result/flnc.fasta.sorted.bam

#04_stat
echo unmap_number=`grep -c 'No paths found for' flnc.fasta.sam.log`
echo unmap_number=`samtools view -f4 flnc.fasta.sorted.bam | wc -l`
python /RSeQC-2.6.6/scripts/bam_stat.py \
-i $workdir/result/flnc.fasta.sorted.bam \
-q 0 > flnc.fasta.sorted.bam_stat.log

```

##### 2) minimap2

```shell
reference= reference/path

#01_index
minimap2 -d hg38.mmi \
$reference/hg38/hg38.fa

#02_minimap
minimap2 -t 10 \
-ax splice \
--MD -uf --secondary=no \
-C5 -O6,24 -B4 hg38.mmi \
$workdir/data/flnc.fasta > $workdir/result/flnc.fasta.sam \
2> flnc.fasta.sam.log

#03_sort
samtools view -bS $workdir/result/flnc.fasta.sam > $workdir/result/flnc.fasta.bam
samtools sort $workdir/result/flnc.fasta.bam > $workdir/result/flnc.fasta.sorted.bam
samtools index $workdir/result/flnc.fasta.sorted.bam

#04_stat
echo unmap_number=`samtools view -f4 flnc.fasta.sorted.bam | wc -l`
python /software/RSeQC-2.6.6/scripts/bam_stat.py \
-i $workdir/result/flnc.fasta.sorted.bam \
-q 0 > flnc.fasta.sorted.bam_stat.log
```

### 基因和转录本鉴定

```shell
#01_tama_collapse
python \
/software/tama/tama_collapse.py \
-b BAM \
-s $workdir/result/flnc.fasta.sorted.bam \
-f $reference/hg38/hg38.fa \
-p isoform2 \
-x no_cap \
-e longest_ends > tama.ollapse.log

```

#### 2) cDNA_Cupcake

```shell
## gmap比对
samtools view $workdir/result/flnc.fasta.sorted.bam \
> $workdir/result/flnc.fasta.sorted.sam
python /software/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py \
 --input $workdir/data/flnc.fasta \
 -s flnc.fasta.sorted.sam \
 --min-coverage 0.99 \
 --min-identity 0.95 \
 --max_5_diff 1000 \
 --max_3_diff 100 \
 --max_fuzzy_junction 10 \
 -o gmap
rm -f $workdir/result/flnc.fasta.sorted.sam

## minimap2比对
samtools view $workdir/result/flnc.fasta.sorted.bam \
> $workdir/result/flnc.fasta.sorted.sam
python /software/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py \
 --input $workdir/data/flnc.fasta \
 -s $workdir/result/flnc.fasta.sorted.sam \
 --min-coverage 0.99 \
 --min-identity 0.95 \
 --max_5_diff 1000 \
 --max_3_diff 100 \
 -o minimap2
rm -f $workdir/result/flnc.fasta.sorted.sam

```



#### 3) SQANTI3(未完)

```shell
export PATH=/local_data1/rna_training/liwei/software/anaconda3/bin:$PATH
export PYTHONPATH=$PYTHONPATH:/local_data1/rna_training/liwei/software/cDNA_Cupcake-9.1.1/sequence/
export PYTHONPATH=$PYTHONPATH:/local_data1/rna_training/liwei/software/cDNA_Cupcake-9.1.1/
source activate SQANTI3.env
python /local_data1/rna_training/liwei/software/SQANTI3/sqanti3_qc.py \
 /local_data1/rna_training/liwei/07.Loci_isoform/00.data/test_chr13_seqs.fasta \
 /local_data1/rna_training/liwei/07.Loci_isoform/00.data/Homo_sapiens.GRCh38.86.chr13.gtf \
 /local_data1/rna_training/liwei/07.Loci_isoform/00.data/Homo_sapiens.GRCh38.dna.chromosome.13.fa \
 --fl_count /local_data1/rna_training/liwei/07.Loci_isoform/00.data/chr13_FL.abundances.txt \
 -o Sample

```

