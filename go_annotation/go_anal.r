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

#加载背景库文件：genes.new.go.annotation 
go_anno <- read.delim('genes.new.go.annotation', header=FALSE, stringsAsFactors =FALSE) 
names(go_anno) <- c('gene_id','ID')
print("背景库加载完成")  

#加载GO注释描述：go_term.list
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

#差异基因与背景库检查 
common_genes <- sum(gene_select %in% go_anno$gene_id)
if (common_genes == 0) {
  stop("差异基因与背景库无交集，请检查ID格式！")
}
print(paste("共有", common_genes, "个差异基因存在于背景库中。"))

#富集分析 
go_rich <- enricher(gene = gene_select,
                    TERM2GENE = go_anno[c('ID','gene_id')],
                    TERM2NAME = go_anno[c('ID','Description')],
                    pvalueCutoff = 0.05,
                    pAdjustMethod = 'BH',
                    qvalueCutoff = 0.2,
                    maxGSSize = 10) 
print("富集分析完成")

dir.create("R_Results")
ENRICH.df = as.data.frame(go_rich)
write.csv(ENRICH.df,"R_Results/go_rich.csv")
message("已保存GO富集文件至: R_Results/go_rich.csv")

if(nrow(go_rich) == 0) {
  print("没有显著富集的GO条目，请调整参数！")
  print(paste("输入基因数:", length(gene_select)))
  print(paste("背景基因数:", length(unique(go_anno$gene_id))))
  stop("富集分析未找到显著结果，请检查输入数据或调整参数。")
}

#类似的可以做另外两个类型的GO富集 
#出图的话，clusterProfiler的可视化也是非常简单 
barplot(go_rich,drop=T)

output_file <- "R_Results/GO_barplot.pdf"  # 可替换为 .png/.tiff/.svg 等
print("成功输出")
# 动态调整图像高度（基于条目数量）
n_cats <- min(20, nrow(go_rich))  # 最多显示20条
img_height <- max(4, n_cats * 0.5) # 基础高度4英寸，每条目增加0.5英寸

# 打开图形设备
pdf(output_file, width = 10, height = img_height)  # PDF格式
# 绘制图形（添加标题）
print(barplot(go_rich, drop = TRUE, 
              title = "GO Enrichment Analysis",
              showCategory = n_cats))

# 关闭设备（保存文件）
dev.off()

message("已保存GO富集条形图至: ", output_file)


#如果要画GO的网络关系图，则需要借助topGO跟Rgraphviz 
#BiocManager::install('topGO') 
#BiocManager::install('Rgraphviz') 
library(topGO) 

write.table(go_rich, 'go_tmp.txt', sep='\t', row.names = FALSE, quote = FALSE) 
tmp <- read.delim('go_tmp.txt') 
tmp <- merge(tmp, go_class[c('ID', 'Ontology')], by = 'ID') 
tmp <- tmp[!is.na(tmp$Ontology), ]  # 删除NA的行
tmp$Ontology <- factor(tmp$Ontology, 
                       levels = c("biological_process", "molecular_function", "cellular_component"),
                       labels = c("BP", "MF", "CC")) %>% as.character()
tmp <- tmp[c(13,1:9)] 
tmp <- tmp[order(tmp$pvalue),] 

generate_GO_plot <- function(ontology, tmp, go_rich_target, 
                             sig_file = "go_rich.significant.txt",
                             width = 10, height = 10, res = 600) {
  
  
  write.table(tmp, sig_file, sep = '\t', row.names = FALSE, quote = FALSE) 
  GO_RICH <- go_rich_target 
  GO_RICH@result <- GO_RICH@result[as.vector(subset(tmp, Ontology == ontology)$ID),] 
  GO_RICH@ontology <- ontology
  
  ###===以下为待添加的if函数===###
  count_entries <- subset(tmp, Ontology == ontology)
  if (nrow(count_entries) == 0) {
    message(paste0("No ", ontology, " ontology terms found. Exiting."))
  }else {
    filename <- paste0("R_Results/GO_Rich_", ontology, ".png")  # 配置图片文件名
    # 绘制GO图并保存
    png(filename, width = width, height = height, units = "in", res = res)
    plotGOgraph(GO_RICH)
    dev.off()
    message("输出图片到：",filename)
    plotGOgraph(GO_RICH)  # 控制台也输出一次
  }
}

for (onto in c("BP", "CC", "MF")) {
  generate_GO_plot( ontology = onto, tmp, go_rich)
}
