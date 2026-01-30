#!/bin/bash
workPATH=[workpath] #工作目录

referencePATH=${workPATH}0_reference/ #参考基因组目录
genome=${referencePATH}GCF_004153795.1_genomic.fasta
GTF=${referencePATH}genomic.gtf
# conda环境激活
condaActPATH=$(whereis activate | awk '{print $2}')
source ${condaActPATH} bio #conda环境激活


# 0. 自动检测可用 CPU 核心数
TOTAL_CORES=$(nproc)
USING_CORES=$((TOTAL_CORES-2))
echo "总核心数: ${TOTAL_CORES}"
# 任务线程分配
    # 如果核心数多，优先分配给并发任务，每个任务给 2-4 线程
    # 核心数大于32，留出 2 个核心给系统，其余按 4 线程/任务切分
    # 如果核心数少，则降低并发，提高单任务线程
if [ ${TOTAL_CORES} -le 4 ]; then
    JOBS=2
    THREADS_PER_JOB=2
elif [ ${TOTAL_CORES} -le 16 ]; then
    JOBS=4
    THREADS_PER_JOB=4
else
    THREADS_PER_JOB=4
    JOBS=$(( (TOTAL_CORES - 2) / THREADS_PER_JOB ))
fi

# 1. fastp
fastpOutPATH=${workPATH}1_fastp_output/ #fastp Output目录
fastpRepoPATH=${workPATH}1_fastp_reports/ #fastp Reports目录
mkdir -p ${workPATH}1_fastp_output
mkdir -p ${workPATH}1_fastp_reports

# 获取所有样本名称
cd ${workPATH}0_originalData
samples=$(ls *.fq.gz | sed 's/_[12]\.fq\.gz//' | sort -u)
echo "找到以下样本:"
echo "$samples"
echo ""

echo "开始处理fastp..."
echo ""
parallel --will-cite -j ${JOBS} '
    sample={}
    echo "处理样本: $sample"
    if [[ -f "${sample}_1.fq.gz" && -f "${sample}_2.fq.gz" ]]; then
        fastp \
            -i "${sample}_1.fq.gz" \
            -I "${sample}_2.fq.gz" \
            -o "${sample}_1.clean.fq.gz" \
            -O "${sample}_2.clean.fq.gz" \
            -h "${sample}.html" \
            -j "${sample}.json" \
            --detect_adapter_for_pe \
            --overrepresentation_analysis \
            --thread ${THREADS_PER_JOB}       
        mv ${sample}*.clean.fq.gz ${fastpOutPATH}
        mv ${sample}.html ${fastpRepoPATH}
        mv ${sample}.json ${fastpRepoPATH}
        echo "完成: $sample"
    else
        echo "错误: 样本 $sample 的配对文件不存在 ----fastp"
    fi
' ::: $samples
echo ""
echo "所有样本处理完成"
echo "结果保存在 fastp_output/ 目录中"
echo "报告保存在 fastp_reports/ 目录中"
cd ${workPATH}

# 2.0 hisat2索引建立
genome_prefix=${referencePATH}GCF_004153795.1_genome_prefix
prefixLogs=${referencePATH}logs/
mkdir -p ${prefixLogs}
hisat2-build ${genome} ${genome_prefix} -p ${TOTAL_CORES} 1>${prefixLogs}hisat2-build.log 2>${prefixLogs}hisat2-build.err

# 2.1 hisat2比对
hisatPATH=${workPATH}hisat2/ #hisat2目录
fastpPATH=${workPATH}1_fastp_output/ #fastp目录
genome_prefix=${referencePATH}prefixOut/GCF_004153795.1_genome_prefix #reference路径
#echo ${genome_prefix} 

mkdir -p ${hisatPATH}
for sample in ${fastpPATH}*_1.clean.fq.gz; do
	sampleName=$(basename ${sample} _1.clean.fq.gz)
    if [ -f ${fastpPATH}${sampleName}_2.clean.fq.gz ]; then		
        hisat2 --new-summary \
        -p 40 \
        -x ${genome_prefix} \
        -1 ${fastpPATH}${sampleName}_1.clean.fq.gz \
        -2 ${fastpPATH}${sampleName}_2.clean.fq.gz \
        -S ${hisatPATH}${sampleName}.sam 2> ${hisatPATH}${sampleName}.err;
	else
		echo -e "错误: 样本 $sample 的R2不存在 ----hisat2";
	fi
done

# 3. 比对结果压缩、排序、构建索引
hisatPATH=${workPATH}2_hisat2/
sortPATH=${workPATH}3_sorted/
mkdir -p ${workPATH}3_sorted/

cd ${hisatPATH}
samples=$(ls *.sam | sed 's/\.sam$//' | tr '\n' ' ')
echo "Samples $samples will be processed"

for file in *.sam; do
    fileName=${file%.sam}
    echo "Processing ${file}..."
    samtools view --threads 10 -b ${file} > ${sortPATH}${file%.sam}.bam   #将SAM格式转换为BAM格式，并删除sam后缀
    echo "Sorting ${fileName}"
    samtools sort -@ 20 -m 2G ${sortPATH}${fileName}.bam > ${sortPATH}${fileName}.sorted.bam   #排序
    echo "Building index for ${fileName}"
    samtools index ${sortPATH}${fileName}.sorted.bam   #创建索引
done

# 4. 计算Count、FPKM、TPM
sortPATH=${workPATH}3_sorted/
feaPATH=${workPATH}4_featureCounts/
sigOutputDir=${feaPATH}singleCount/
mkdir -p ${feaPATH}
mkdir -p ${sigOutputDir}

# 多线程执行
    # cd ${feaPATH}
    # parallel -j ${JOBS} \
    #   "Rscript ${workPATH}4_featureCounts.R {} ${referencePATH}genomic.gtf ${USING_CORES} ${feaPATH}singleCount/{/.}" \
    #   ::: ${sortPATH}*.sorted.bam
    # cd ${workPATH}

cd ${feaPATH}
for file in ${sortPATH}*.sorted.bam; do
    sampleName=$(basename ${file} .sorted.bam)
    Rscript ${workPATH}4_featureCounts.R "$file" ${GTF} ${USING_CORES} "${sigOutputDir}${sampleName}"; 
done
echo "计算结果保存在singleCount"
cd ${workPATH}

# 合并Count
perl -MFile::Basename -lne '
    $filename = basename($ARGV, ".count");
    print "$filename\t$_"
' ${feaPATH}singleCount/*.count > ${feaPATH}merge.count
echo "保存了合并文件：merge.count"
#生成矩阵
awk -F"\t" '{print $1"\t"$2"\t"$3}' ${feaPATH}merge.count > ${feaPATH}gene.FPKMcount
awk -F"\t" '{print $1"\t"$2"\t"$5}' ${feaPATH}merge.count > ${feaPATH}gene.tpm
cd ${feaPATH}
# R和shell分离
    # Rscript ${workPATH}matrixCountAndTPM.R
#嵌入R
R --vanilla --slave <<EOF
    library(reshape2)
    # 读取 mergeCount
    geneCount=read.table('gene.FPKMcount',header=FALSE,sep="\t")
    geneTPM=read.csv('gene.tpm',header=F,sep="\t")
    # 设置列名
    colnames(geneCount) = c('sample', 'gene', 'counts')
    colnames(geneTPM)=c('sample','gene','tpm')
    # 设置格式
    transCount = dcast(geneCount, formula = gene ~ sample, value.var = "counts", fill = 0)
    transTPM = dcast(geneTPM,formula=gene~sample)
    # 输出
    write.table(transCount,file="matrixFPKMcount.txt", sep="\t", quote=FALSE, row.names=FALSE)
    write.table(transTPM,file="matrixTpm.txt",sep="\t",quote=FALSE,row.names=FALSE)
EOF

echo "保存了：matrixFPKMcount.txt、matrixTpm.txt"
cd ${workPATH}
