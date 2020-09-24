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
Trinity --seqType fq --max_memory 100G --no_salmon \                                            ##由于Trinity版本问题，有时salmon与之不兼容，所以添加--no_salmon；
  --left $workdir/DAYCA_result/clean/DAYCA_1_clean.fq.gz \
  --right $workdir/DAYCA_result/clean/DAYCA_2_clean.fq.gz \
  --CPU 20 \
  --output $workdir/DAYCA_result/trinity2_out > $workdir/workplace/log/DAYCA_trinity.log 2>&1
```

## 第四步：blastn比对到病毒数据库
