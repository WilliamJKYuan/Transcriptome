#install.packages("BiocManager")

#检查运行环境
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
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("clusterProfiler")


#加载clusterProfiler 
library(clusterProfiler)
  
#加载背景库文件，即我们的第二个文件 
go_anno <- read.delim('genes.new.go.annotation', header=FALSE, stringsAsFactors =FALSE) 
names(go_anno) <- c('gene_id','ID')
go_anno$gene_id <- sub("\\.p\\d+$", "", go_anno$gene_id)  # 移除.pX后缀
print("背景库加载完成")  

#加载GO注释描述，即我们的第三个文件 
go_class <- read.delim('go_term.list', header=FALSE, stringsAsFactors =FALSE) 
names(go_class) <- c('ID','Description','Ontology')
print("GO注释描述加载完成")  

#合并背景与GO文件 
go_anno <-merge(go_anno, go_class, by = 'ID', all.x = TRUE)  
print("合并背景与GO文件完成")  

#差异基因导入 
gene_list <- read.delim('gene.list',stringsAsFactors = FALSE) 
names(gene_list) <- c('gene_id') 
gene_select <- gene_list$gene_id
print("差异基因导入完成")

#富集分析 
go_rich <- enricher(gene = gene_select,
                  TERM2GENE = go_anno[c('ID','gene_id')],
                  TERM2NAME = go_anno[c('ID','Description')],
                  pvalueCutoff = 0.5,
                  pAdjustMethod = 'BH',
                  qvalueCutoff = 2.5,
                  maxGSSize = 2500) 
print("富集分析完成")

ENRICH.df = as.data.frame(go_rich)
write.csv(ENRICH.df,"go_rich.csv")
message("已保存GO富集文件至: go_rich.csv")

if(nrow(go_rich) == 0) {
  print("没有显著富集的GO条目，请调整参数！")
  print(paste("输入基因数:", length(gene_select)))
  print(paste("背景基因数:", length(unique(go_anno$gene_id))))
  stop("富集分析未找到显著结果，请检查输入数据或调整参数。")
}

#使用clusterProfiler进行可视化
barplot(go_rich,drop=T)
output_file <- "GO_barplot.pdf"  # 可替换为 .png/.tiff/.svg 等
print("成功输出")
# 动态调整图像高度（基于条目数量）
n_cats <- min(20, nrow(go_rich))  # 最多显示20条
img_height <- max(4, n_cats * 0.5) # 基础高度4英寸，每条目增加0.5英寸
# 打开图形设备
pdf(output_file, width = 10, height = img_height)  # PDF格式
# 绘制图形
print(barplot(go_rich, drop = TRUE, 
              title = "GO Enrichment Analysis",
              showCategory = n_cats))
dev.off()
message("已保存GO富集条形图至: ", output_file)


#如果要画GO的网络关系图，则需要借助topGO跟Rgraphviz 
#BiocManager::install('topGO') 
#BiocManager::install('Rgraphviz') 
library(topGO) 

write.table(go_rich, 'go_tmp.txt', sep='\t', row.names = FALSE, quote = FALSE) 
tmp <- read.delim('go_tmp.txt') 
tmp <- merge(tmp, go_class[c('ID', 'Ontology')], by = 'ID') 
tmp <- tmp[c(10,1:9)] 
tmp <- tmp[order(tmp$pvalue),] 

write.table(tmp, 'go_rich.significant.txt', sep = '\t', row.names = FALSE, quote = FALSE) 
go_rich_BP <- go_rich 
go_rich_BP@result <- go_rich_BP@result[as.vector(subset(tmp,  Ontology == 'biological process')$ID),] 
go_rich_BP@ontology <- 'BP' 
plotGOgraph(go_rich_BP)
