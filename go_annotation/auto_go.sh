source /home/yuan/enter/bin/activate trans
#!/bin/bash
mkdir Results

#保留最长转录本
TransDecoder.LongOrfs -t NoEnter.Unigene.fa
mv ./NoEnter.Unigene.fa.transdecoder_dir/longest_orfs.cds ./Results

#Swiss建库
diamond makedb --in swissprot -d swissprot
echo "Build Success"

#使用diamond对核酸序列进行比对
echo "核酸序列比对"
diamond blastx -d swissprot -q ./Results/longest_orfs.cds -k 1 -e 0.00001 -o ./Results/genes.swiss_dia_matches.m8
echo "比对完成"

转换成GO_ID
#echo "开始转换GO_ID"
python ./py/get_swiss_go.py idmapping.tb.gz > swiss_go.list
echo "转换完成"

#提取diamond比对结果
echo "提取diamond比对结果"
python ./py/get_trinity_swiss_id.py ./Results/genes.swiss_dia_matches.m8 > ./Results/trinity_swiss.id
echo "提取完成"

source /home/yuan/enter/bin/deactivate
source /home/yuan/enter/bin/deactivate
echo "conda环境已关闭"

#GO注释
echo "开始GO注释"
python -u ./py/get_go_annotation.py swiss_go.list ./Results/trinity_swiss.id | tee

#go-basic_ID解析
python ./py/get_go_term.py go-basic.obo

#修改格式
python ./py/split_with_one_go.py ./Results/genes.go.annotation ./Results/genes.new.go.annotation

#R语言GO分析
Rscript go_anal.r

echo "Cleaning Up..."
mv ./NoEnter.Unigene.fa.transdecoder_dir ./Results
mv go_rich.significant.txt ./R_Results
mv go_tmp.txt ./R_Results
echo "All Done ( •̀ ω •́ )y"
