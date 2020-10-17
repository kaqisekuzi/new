## 第一步：fastqc检测raw data的质量
```
fastqc -o /share/home/stu_marui/insect/result/clean_data -f fastq ~/insect/data/ZLA_1.fq.gz
fastqc -o /share/home/stu_marui/insect/result/clean_data -f fastq ~/insect/data/ZLA_2.fq.gz
```

## 第二步：fastp 进行质量控制
```
workdir=/share/home/stu_marui/insect

fastp -i $workdir/data/LYJA_1.fq.gz \
      -I $workdir/data/LYJA_2.fq.gz \
      -o $workdir/LYJA_result/clean/LYJA_1_clean.fq.gz \
      -O $workdir/LYJA_result/clean/LYJA_2_clean.fq.gz \
      -W 5 -M 20 -5 -3 -l 50 -w 10 -j $workdir/workplace/stat/LYJA_QC_5.json > $workdir/workplace/log/LYJA_QC.log 2>&1
```

## 第三步：Trinity 组装
```
workdir=/share/home/stu_marui/insect

## 由于Trinity版本问题，有时salmon与之不兼容，所以添加--no_salmon；
Trinity --seqType fq --max_memory 100G --no_salmon \                                            
  --left $workdir/DAYCA_result/clean/DAYCA_1_clean.fq.gz \
  --right $workdir/DAYCA_result/clean/DAYCA_2_clean.fq.gz \
  --CPU 20 \
  --output $workdir/DAYCA_result/trinity2_out > $workdir/workplace/log/DAYCA_trinity.log 2>&1
```

## 第四步：blastn比对到病毒数据库
# virus是提前建好的库
```
workdir=/share/home/stu_marui/insect
virus=/share/home/stu_marui/Pst_RNAseq/RNA_seq_Cicadellidae/workspace/virus

blastn -db $virus/virus_blastdb \
       -max_target_seqs 1 \
       -query $workdir/DAYCA_result/trinity2_out/Trinity.fasta \
       -outfmt 6 \
       -out $workdir/DAYCA_result/DAYCA_trinity2virus.fmt6
```

## 第五步：make annotation file
```
workdir=/share/home/stu_marui/insect
virus=/share/home/stu_marui/Pst_RNAseq/RNA_seq_Cicadellidae/workspace/virus

grep "^>" $virus/virus.fasta |
    awk -v FS="|" -v OFS="\t" '{a=substr($1,2);sub(" *$", "", a);print a, $2}' \
    > $workdir/virus_genome_annotation.txt
```

## 第六步：tidy blast result
```
workdir=/share/home/stu_marui/insect

python $workdir/workplace/tidy_blast_result.py $workdir/CXA_result/CXA_trinity2virus.fmt6 \
       $workdir/virus_genome_annotation.txt \
       -c $workdir/CXA_result/blast2virus.count.csv $workdir/CXA_result/blast2virus.csv
```
#### tidy_blast_result.py
#####################################################################################################################
```

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# Author       : _
# CreateTime   : 2020-09-11 15:02:20
# LastEditTime : 2020-09-11 16:58:02
# LastEditors  : Keliang_Lyu
# Description  : Tidy blast result and add annotation.
# ====================================
import argparse
import pandas as pd


def load_fmt6(fmt6path):
    fmt6_header = [
        'Queryid', 'Targetid', 'Identity', 'Alignment_length', 'Mismatch',
        'Gap_open', 'Q_start', 'Q_end', 'S_start', 'S_end', 'Evalue', 'Score'
    ]
    blastdf = pd.read_csv(fmt6path, sep='\t', names=fmt6_header)
    return blastdf


def anno_blast(blastdf, annodf):
    resultheader = [
        'Queryid', 'Targetid', 'Annotation', 'Identity', 'Alignment_length',
        'Mismatch', 'Gap_open', 'Q_start', 'Q_end', 'S_start', 'S_end'
    ]
    resultdf = blastdf.merge(
        annodf, how='left', on='Targetid'
    )
    resultdf = resultdf.loc[:, resultheader]
    return resultdf


def tidy_blast_result(fmt6path, annopath):
    blastdf = load_fmt6(fmt6path)
    annodf = pd.read_csv(annopath, sep='\t', names=('Targetid', 'Annotation'))
    resultdf = anno_blast(blastdf, annodf)
    resultdf.sort_values('Alignment_length', ascending=False, inplace=True)
    return resultdf


def count_anno(resultdf):
    countdf = pd.DataFrame(resultdf['Targetid'].value_counts())
    countdf.reset_index(inplace=True)
    countdf.columns = ['Targetid', 'Count']
    annodf = resultdf.loc[:, ['Targetid', 'Annotation']]
    annodf = annodf.drop_duplicates()
    countdf = countdf.merge(annodf, how='left', on='Targetid')
    countdf = countdf.loc[:, ['Targetid', 'Annotation', 'Count']]
    countdf.sort_values('Count', ascending=False, inplace=True)
    return countdf


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Blast result file with format 6.')
    parser.add_argument('anno', help='Annotation file.')
    parser.add_argument('output', help='Output file.')
    parser.add_argument('-c', '--count', help='Output Count file.')
    args = parser.parse_args()

    resultdf = tidy_blast_result(args.input, args.anno)
    resultdf.to_csv(args.output, index=False)

    # Count annotation
    if args.count:
        countdf = count_anno(resultdf)
        countdf.to_csv(args.count, index=False)
```
#######################################################################################################################




