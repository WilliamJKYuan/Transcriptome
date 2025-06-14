### 加载R包
library(GenomicFeatures)
library(DESeq2)
library(dplyr)

### 设置输入和输出文件路径
setwd("{Project_Path}")
file_in <- "All_gene_expression_countData.txt"        # 输入的基因计数数据文件
file_design <- "All_gene_expression_groupInfo.txt"    # 数据分组文件
file_compare <- "compare.txt"                         # 比较组信息文件
file_deg_num = paste("data//", "DE_", file_in, sep="")  ## 所有比较中的差异表达基因数
file_final_csv  = paste("data//", "DE_", file_in, "_Final_Out.csv", sep="")  ## 所有基因的最终输出
file_final_genelist = paste("data//", "DEG_geneid_allcomapre.txt", sep="")   ## 所有差异基因ID列表

#创建输出目录
# 创建 'data' 文件夹，如果它不存在的话
if (!dir.exists("data")) {
  dir.create("data")
}

# 创建 'data/DEGsid' 文件夹，如果它不存在的话
if (!dir.exists("data/DEGsid")) {
  dir.create("data/DEGsid")
}

#######获取不同比较中的差异表达基因 (DEGs)

# 读取计数数据文件
data_in = read.table(file_in, head=TRUE, row.names=1, check.names=FALSE)

# 读取比较信息文件
mycompare = read.table(file_compare, head=TRUE)

# 读取实验设计文件
mydesign = read.table(file_design, head=TRUE)

## 过滤掉总和为0的基因计数行
countData = as.data.frame(data_in)
dim(countData)  # 输出计数数据的维度
mycounts_filter <- countData[rowSums(countData) != 0,]
dim(mycounts_filter)  # 输出过滤后的计数数据的维度

## 检查比较组信息和实验设计
head(mycompare)
head(mydesign)

# 获取比较组数目
total_num = dim(mycompare)[1]  
tracking = 0
gene_num_out = c()
pvalue_cut = 0.01  ## P值阈值
condition_name = c()

# 对每一组比较进行循环处理
for (index_num in c(1:total_num)) {
  
  tracking = tracking + 1
  
  test_name = as.character(mycompare[index_num, 1])  # 获取当前比较的名称
  group1 = as.character(mycompare[index_num, 2])     # 获取组1
  group2 = as.character(mycompare[index_num, 3])     # 获取组2
  
  # 根据组名从实验设计中提取样本信息
  sample1 = mydesign[mydesign$Group == group1,]
  sample2 = mydesign[mydesign$Group == group2,]
  allsample = rbind(sample1, sample2)
  
  # 获取样本ID和组别信息
  counts_sample = as.character(allsample$counts_id)
  groupreal = as.character(allsample$Group)
  
  # 提取对应样本的基因计数数据
  countData = mycounts_filter[, counts_sample]
  colData = as.data.frame(cbind(counts_sample, groupreal))
  names(colData) = c("sample", "condition")
  
  # 创建DESeq2数据集
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
  
  # 运行DESeq分析
  dds <- DESeq(dds)
  
  # 获取差异分析的结果
  resSFtreatment <- results(dds, cooksCutoff = FALSE, contrast = c("condition", group2, group1))
  
  # 将结果转为数据框
  out_test = as.data.frame(resSFtreatment)
  
  # 提取log2FoldChange、p-value和adjusted p-value
  final_each = cbind(out_test$log2FoldChange, out_test$pvalue, out_test$padj)
  rownames(final_each) = resSFtreatment@rownames
  
  # 给列命名
  names = c('(logFC)', '(pvalue)', '(Qvalue)')
  final_name = paste(test_name, '_', names, sep="")
  colnames(final_each) = final_name
  
  # 如果是第一个比较组，初始化最终表格
  if (tracking == 1) {
    final_table = final_each
  } else {
    final_table = cbind(final_table, final_each)
  }
  
  # 筛选符合条件的基因：P值小于0.05，log2FoldChange绝对值大于等于1，padj小于0.05
  gene_sel = out_test[((!is.na(out_test$pvalue)) & (!is.na(out_test$log2FoldChange))) & out_test$pvalue < 0.05 & abs(out_test$log2FoldChange) >= 1 & out_test$padj < 0.05,]
  
  gene_sel <- na.omit(gene_sel)  ## 删除包含NA的行
  
  # 分别提取上调和下调基因
  gene_sel_up = gene_sel[gene_sel$log2FoldChange > 0,]
  gene_sel_do = gene_sel[gene_sel$log2FoldChange < 0,]
  
  # 输出上调和下调基因的文件
  file_out_up = paste("data//DEGsid//", "UP_", test_name, ".csv", sep="")
  file_out_do = paste("data//DEGsid//", "DOWN_", test_name, ".csv", sep="")
  
  # 提取上调和下调基因的基因ID列表
  gene_list_up = rownames(gene_sel_up)
  gene_list_do = rownames(gene_sel_do)
  
  # 将上调和下调基因合并成列表
  all_ub_down = list(gene_list_up, gene_list_do)
  
  nameup <- paste(test_name, "_up", sep="")
  namedown <- paste(test_name, "_down", sep="")
  names(all_ub_down) = c(nameup, namedown)
  
  # 如果是第一个比较组，初始化最终的基因列表
  if (tracking == 1) {
    final_genelist = all_ub_down
  } else {
    final_genelist = c(final_genelist, all_ub_down)
  }
  
  # 计算上调和下调基因的数量
  gene_num_up = length(gene_list_up)
  gene_num_do = length(gene_list_do)
  
  # 将上调和下调基因写入文件
  write.table(gene_sel_up, file = file_out_up, row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
  write.table(gene_sel_do, file = file_out_do, row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)
  
  # 更新条件名称和基因数量
  condition_name = c(condition_name, paste("UP_", test_name, sep=""), paste("DO_", test_name, sep=""))
  gene_num_out = c(gene_num_out, gene_num_up, gene_num_do)
}

# 将最终的差异基因列表合并
final_DEGs_list <- do.call(cbind, lapply(lapply(final_genelist, unlist), `length<-`, max(lengths(final_genelist))))

# 输出最终的差异基因列表
final_DEGs_list

# 合并条件名称和差异基因数量
out_final2 = cbind(condition_name, gene_num_out)
colnames(out_final2) = c("Tests", "DEG number")

# 将结果写入输出文件
write.table(final_DEGs_list, file = file_final_genelist, row.names = FALSE, sep = ",", na = "", quote = FALSE)   ## 所有差异基因列表

write.table(out_final2, file = file_deg_num, row.names = FALSE, sep = ",", quote = FALSE)  ## 各比较中的差异基因数量

write.csv(final_table, file = file_final_csv, row.names = TRUE, quote = TRUE)  ## 所有基因的最终输出
