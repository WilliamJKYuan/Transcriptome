####### 检查r包 #######
#pkg install if missing
piif <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(
      {
        install.packages(pkg)
        # 安装后二次验证
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop("包安装失败: ", pkg)
        }
      },
      error = function(e) {
        stop("安装过程中出错: ", e$message)
      }
    )
  }
}
piif("GenomicFeatures")
piif("DESeq2")
piif("dplyr")
piif("foreach")
piif("doParallel")
piif("parallel")

### 加载R包
library(GenomicFeatures)
library(DESeq2)
library(dplyr)
library(foreach)
library(doParallel)
library(parallel)

### 设置输入和输出文件路径
setwd("{Project_Path}")
file_in <- "All_gene_expression_countData.txt"        # 输入的基因计数数据文件
file_design <- "All_gene_expression_groupInfo.txt"    # 数据分组文件
file_compare <- "compare.txt"                         # 比较组信息文件
file_deg_num = paste("data//", "DE_", file_in, sep="")  ## 所有比较中的差异表达基因数
file_final_csv  = paste("data//", "DE_", file_in, "_Final_Out.csv", sep="")  ## 所有基因的最终输出
file_final_genelist = paste("data//", "DEG_geneid_allcomapre.txt", sep="")   ## 所有差异基因ID列表

# 创建输出目录
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("data/DEGsid")) dir.create("data/DEGsid")

####### 数据预处理 ####### 
# 读取数据
## 读取计数数据文件
data_in <- read.table(file_in, header=TRUE, row.names=1, check.names=FALSE)
## 读取比较信息文件
mycompare <- read.table(file_compare, header=TRUE)
## 读取实验设计文件
mydesign <- read.table(file_design, header=TRUE)

# 过滤掉总和为0的基因计数行
countData <- as.data.frame(data_in)
dim(countData)  # 输出计数数据的维度
mycounts_filter <- countData[rowSums(countData) != 0, ]
dim(mycounts_filter)  # 输出过滤后的计数数据的维度

#结果检查
## 检查比较组信息和实验设计
head(mycompare)
head(mydesign)
## 获取比较组数目
total_num = dim(mycompare)[1]  
tracking = 0
gene_num_out = c()
pvalue_cut = 0.01  ## P值阈值
condition_name = c()

####### 并行计算设置 #######
num_cores <- detectCores() - 1  # 留一个核心给系统
cl <- makeCluster(num_cores)
registerDoParallel(cl)
print(paste("可用核心",num_cores))

####### 定义并行处理函数 #######
process_contrast <- function(index_num, mycompare, mydesign, mycounts_filter) {
  test_name <- as.character(mycompare[index_num, 1])    # 获取当前比较的名称
  group1 <- as.character(mycompare[index_num, 2])       # 获取组1
  group2 <- as.character(mycompare[index_num, 3])       # 获取组2
  print(paste("开始比较",test_name))
  
  # 获取样本ID和组别信息
  samples <- rbind(
    mydesign[mydesign$Group == group1, ],
    mydesign[mydesign$Group == group2, ]
  )
  
  # 提取对应样本的基因计数数据
  countData <- mycounts_filter[, samples$counts_id]
  colData <- data.frame(
    sample = samples$counts_id,
    condition = samples$Group
  )
  
  # DESeq2分析流程
  ## 创建DESeq2数据集
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
  ## 运行DESeq分析
  dds <- DESeq(dds)
  ## 获取差异分析的结果
  res <- results(dds, contrast = c("condition", group2, group1), cooksCutoff = FALSE)
  
  # 提取结果
  final_each <- data.frame(
    logFC = res$log2FoldChange,
    pvalue = res$pvalue,
    padj = res$padj
  )
  colnames(final_each) <- paste(test_name, c("logFC", "pvalue", "padj"), sep = "_")
  rownames(final_each) <- rownames(res)
  
  # 筛选差异基因
  deg <- subset(as.data.frame(res), 
                pvalue < 0.05 & 
                  abs(log2FoldChange) >= 1 & 
                  padj < 0.05 & 
                  !is.na(pvalue))
  
  # 输出差异基因文件
  up_genes <- rownames(subset(deg, log2FoldChange > 0))
  down_genes <- rownames(subset(deg, log2FoldChange < 0))
  
  write.csv(deg[deg$log2FoldChange > 0, ], paste0("data/DEGsid/UP_", test_name, ".csv"))
  write.csv(deg[deg$log2FoldChange < 0, ], paste0("data/DEGsid/DOWN_", test_name, ".csv"))
  
  # 返回结果组件
  list(
    final_each = final_each,
    deg_counts = c(length(up_genes), length(down_genes)),
    deg_names = c(paste0("UP_", test_name), paste0("DOWN_", test_name)),
    deg_lists = list(up = up_genes, down = down_genes)
  )
}

####### 并行执行主流程 #######
print("处理中")
total_num <- nrow(mycompare)
results <- foreach(i = 1:total_num, 
                   .packages = c("DESeq2"),
                   .export = c("mycompare", "mydesign", "mycounts_filter")) %dopar% {
                     process_contrast(i, mycompare, mydesign, mycounts_filter)
                   }

stopCluster(cl)  # 关闭并行集群

####### 结果整合 #######
# 合并表达量表
final_table <- do.call(cbind, lapply(results, function(x) x$final_each))

# 合并差异基因统计
gene_num_out <- unlist(lapply(results, function(x) x$deg_counts))
condition_name <- unlist(lapply(results, function(x) x$deg_names))

# 合并基因列表
final_genelist <- do.call(c, lapply(results, function(x) x$deg_lists))
final_DEGs_list <- do.call(cbind, lapply(final_genelist, function(x) `length<-`(x, max(lengths(final_genelist)))))

####### 文件输出 #######
# 输出差异基因统计
write.csv(data.frame(Tests = condition_name, DEG_number = gene_num_out), 
          file_deg_num, 
          row.names = FALSE)

# 输出全量结果
write.csv(final_table, file_final_csv)

# 输出基因列表
write.csv(final_DEGs_list, file_final_genelist, na = "")
