#install.packages("dplyr")
#install.packages("tibble")
#install.packages("ggplot2")
#install.packages("xlsx")
#install.packages("Rmisc") 

knitr::opts_chunk$set(warning = F, message = F)
library(dplyr)
library(tibble)
library(ggplot2)
library(xlsx)
library(Rmisc)
#加载R包

get_qPCR <- function(dataset=dat,
                     ref_gene="Actin",
                     control_group="R",
                     grp=c("F1", "M1", "M6")){
  
  # dataset=dat                   # 初始数据
  # ref_gene="GAPDH"              # 内参基因名字
  # control_group="6H NC"         # 对照组
  # grp=c("Group1", "Gropu2")                # 实验组排序

#参数设置
 # 参数检查
  if(!any(is.element(colnames(dataset), c("Sample_Name", "Target_Name", "CT")))){
    stop("Check the sheet's colnames")
  }
  sampleid <- c("Sample_Name", "Target_Name", "CT")
  dat <- dataset %>% select(all_of(sampleid))
  
  # Step1：计算对照组和处理组的内参基因平均值
  
  # Step1.1: 筛选内参基因时同时限制处理组 
  dat_ref_gene <- dat %>% 
    filter(Target_Name == ref_gene, 
           Sample_Name %in% c(control_group, grp))
  
  # Step1.2：计算平均值
  ref_gene_mean <- dat_ref_gene %>% 
    group_by(Sample_Name) %>%
    dplyr::summarise(CT_ref_mean = mean(CT, na.rm=TRUE))
  
  # Step2: 计算对照组和处理组待检测目的基因减去对应分组的内参基因的平均Ct值
  # Step2.1: 处理目的基因数据时限制组别
  dat_gene <- dat %>% 
    filter(Target_Name != ref_gene,
           Sample_Name %in% c(control_group, grp))
  
  # Step2.2: 合并内参均值并计算ΔCt
  dat_gene_merge <- dat_gene %>% 
    inner_join(ref_gene_mean, by="Sample_Name") %>%
    mutate(CT_delta = CT - CT_ref_mean) 
  
  # Step2.3：分离对照组与处理组
  dat_control <- dat_gene_merge %>% 
    filter(Sample_Name == control_group) %>%
    group_by(Sample_Name, Target_Name) %>%
    dplyr::summarise(Delta_CT_control_mean=mean(CT_delta, na.rm=TRUE)) %>% 
    dplyr::rename(Sample_Name_control=Sample_Name)
  
  dat_treat <- dat_gene_merge %>% 
    filter(Sample_Name != control_group) %>%
    dplyr::rename(Sample_Name_treat=Sample_Name)
  
  # Step3: 计算对照组检测基因的平均Δ值，基于对照组检测基因的平均Δ值，计算实验组的2-ΔΔCt值
  dat_double_delta <- inner_join(dat_treat, dat_control, by="Target_Name") %>%
    mutate(CT_delta_delta = CT_delta - Delta_CT_control_mean,
           qPCR = 2^(-CT_delta_delta))
  
  # Step4：绘图数据准备
  dat_plot <- dat_double_delta %>% 
    dplyr::rename(Sample_Name=Sample_Name_treat) %>%
    filter(Sample_Name %in% grp) %>% # 确保仅保留指定处理组
    dplyr::select(Sample_Name, Target_Name, qPCR)
  
  # 统计时自动排除不在grp中的组
  dat_plot_bar <- Rmisc::summarySE(dat_plot, measurevar="qPCR", 
                                  groupvars=c("Sample_Name", "Target_Name")) %>%
    mutate(Sample_Name = factor(Sample_Name, levels=grp), # 设置因子水平
           Target_Name = factor(Target_Name)) %>%
    filter(!is.na(Sample_Name)) # 排除因错误grp输入导致的NA

  # Step5: 条形图或相关性散点图可视化
  dat_plot <- dat_double_delta %>% 
    dplyr::rename(Sample_Name=Sample_Name_treat) %>%
    dplyr::select(Sample_Name, Target_Name, qPCR) 
  dat_plot_bar <- Rmisc::summarySE(dat_plot, measurevar = "qPCR", 
                                   groupvars = c("Sample_Name", "Target_Name")) %>%
    mutate(Sample_Name=factor(Sample_Name, levels = grp),
           Target_Name=factor(Target_Name)) %>% 
    group_by(Sample_Name, Target_Name) %>%
    mutate(ylimit=(qPCR+sd)) %>%
    ungroup()
  
  dat_plot_bar_ymax <- dat_plot_bar %>% 
    group_by(Target_Name) %>% 
    summarise_at(vars(ylimit), max)
  
  # dat_plot_range <- dat_plot %>% group_by(Sample_Name, Target_Name) %>%
  #   summarise(ymin=min(qPCR), ymax=max(qPCR))
  # setting y axis scale
  y_group <- c()
  y_scale <- c()
  for(i in 1:nrow(dat_plot_bar_ymax)){
    y_group <- c(y_group, rep(as.character(dat_plot_bar_ymax$Target_Name[i]), 2))
    y_scale <- c(y_scale, c(0, ceiling(dat_plot_bar_ymax$ylimit[i])))
  }
  blank_data <- data.frame(Target_Name = y_group, 
                           Sample_Name = 1, 
                           qPCR = y_scale)
  
  # step6: visualization
  pl <- ggplot(dat_plot_bar, aes(x=Sample_Name, weight=qPCR))+
    geom_hline(aes(yintercept = qPCR), color = "gray")+
    geom_bar(color = "black", width = .4, position = "dodge")+
    geom_errorbar(aes(ymin = qPCR, ymax = qPCR + se), 
                  width = 0.25, size = 0.5, position = position_dodge(0.7))+
    labs(x="", y=expression(paste(log[2], " fold change in expression")))+ 
    geom_blank(data = blank_data, aes(x = Sample_Name, y = qPCR))+
    expand_limits(y = 0)+
    scale_y_continuous(expand = c(0, 0))+
    facet_wrap(. ~ Target_Name, scales = "free")+
    theme_bw()+
    theme(axis.title = element_text(face = "bold", color = "black", size = 14),
          axis.text = element_text(color = "black", size = 10),
          axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
          text = element_text(size = 10, color = "black", family="serif"),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.height = unit(0.6, "cm"),
          legend.text = element_text(face = "bold", color = "black", size = 10),
          strip.text = element_text(face = "bold", size = 14))
  res <- list(dat=dat_double_delta, plot=pl)
  return(res)  
}

dat <- read.xlsx("0530_ROUND3.xlsx", sheetIndex = 1)
head(dat)

#数据输入
qPCR_res <- get_qPCR()

DT::datatable(qPCR_res$dat)
#计算结果

qPCR_res$plot#可视化
