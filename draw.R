

#相关R包载入：
library(dplyr)
library(ggplot2)
library(ggrepel)
#本地测试数据读入：
df <- rio::import("data2.txt")

colnames(df)[9] <- "qvalue"

# [1] "seqnames"      "start"         "end"           "width"         "strand"        "peak_name"     "summit"        "foldchange"    "log10(qvalue)"
# [10] "annotation"    "SYMBOL"

#阈值确定：
cutoff_qvalue = 6
log2FC = 5
#根据阈值添加上下调分组标签：
df$group <- case_when(
  df$foldchange > log2FC & df$qvalue > cutoff_qvalue ~ "up",
  TRUE ~ 'none'
)
head(df)


#转换为因子，指定绘图顺序；
#df$'qvalue' <- df$qvalue #新增-log10p列
df$group <- factor(df$group, levels = c("up","none"))

#自定义颜色：
mycol <- c("#EB4232","#d8d8d8")
#自定义主题：
mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5))


#ggplot2绘制火山图：
p <- ggplot(data = df,
            aes(x = foldchange,
                y = qvalue,
                color = group)) + #建立映射
  geom_point(size = 2.2) + #绘制散点
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + #自定义散点颜色
  scale_x_continuous(limits = c(0, 140),
                     breaks = c(5, seq(0, 140, by = 20)), # 在5处添加横坐标标签
                     labels = c("5", as.character(seq(0, 140, by = 20)))) +
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 2200),
                     breaks = c(6, seq(0, 2200, by = 200)), # 在6处添加纵坐标标签
                     labels = c("6", as.character(seq(0, 2200, by = 200)))) +
  geom_hline(yintercept = c(cutoff_qvalue),size = 0.7,color = "black",lty = "dashed") + #水平阈值线
  geom_vline(xintercept = c(log2FC),size = 0.7,color = "black",lty = "dashed") + #垂直阈值线
  mytheme
p


#分别筛选上下调中显著性top20，作为本次测试需要添加的目标标签（共40个标签）：
up <- filter(df, group == 'up') %>% distinct(SYMBOL, .keep_all = T) %>%
  top_n(20, qvalue)

up <- filter(df,SYMBOL %in% c("dnd1","nanos3","ddx4","kop","ca15b","gra","rgs14a","tdrd7a","celf1","dazl","hook2"))

#两个函数都可用于添加标签：
##geom_text_repel():常规样式，默认参数配置
p1 <- p +
  geom_text_repel(data = up,
                  aes(x = foldchange, y = qvalue, label = SYMBOL))
p1

##geom_label_repel():带边框样式，默认参数配置
p2 <- p +
  geom_label_repel(data = up,
                   aes(x = foldchange, y = qvalue, label = SYMBOL))
p2

#重点参数学习—以常规样式为例：
p3 <- p +
  geom_point(data = up,
             aes(x = foldchange, y = qvalue),
             color = '#EB4232', size = 4.5, alpha = 0.2) + #图中强调显示一下需要标注的散点
  geom_text_repel(data = up,
                  aes(x = foldchange, y = qvalue, label = SYMBOL),
                  seed = 233, #可设置随机数种子
                  size = 3.5, #字体大小
                  color = 'black', #字体颜色
                  min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
                  force = 2, #重叠标签间的排斥力
                  force_pull = 2, #标签和数据点间的吸引力
                  box.padding = 0.4, #标签周边填充量，默认单位为行
                  max.overlaps = Inf #排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
  )
p3

#不添加指引线段效果：
p +
  geom_point(data = up,
             aes(x = foldchange, y = qvalue),
             color = '#EB4232', size = 4.5, alpha = 0.2) +
  geom_text_repel(data = up,
                  aes(x = foldchange, y = qvalue, label = SYMBOL),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = Inf, #改为Inf去掉指引线段
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.4,
                  max.overlaps = Inf)



#个性化调整：对齐目标标签
p5 <- p +
  geom_point(data = up,
             aes(x = foldchange, y = qvalue),
             color = '#EB4232', size = 4.5, alpha = 0.2) +
  geom_text_repel(data = up,
                  aes(x = foldchange, y = qvalue, label = SYMBOL),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.5, #线段不透明度
                  nudge_x = 140 - up$foldchange, #标签x轴起始位置调整
                  direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 0 #对齐标签：0右对齐，1左对齐，0.5居中
  ) + 
  labs(x = "foldchange", y = "-log10(qvalue)")
p5




# GO and KEGG enrichment analysis
setwd("/home/huangqllab/su20240422/6.callpeak")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(org.Mm.eg.db)

txdb=TxDb.Mmusculus.UCSC.mm10.knownGene

peak <- readPeakFile('only_marf1_peak.bed', as="GRanges")

peakAnno <- annotatePeak(peak,
                         tssRegion = c(-500, 500),
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")

annotation <- as.data.frame(peakAnno)
out_name="/home/huangqllab/su20240422/7.annotation/only_marf1.txt"
rio::export(annotation,file=out_name,col.names=T)




library(clusterProfiler)
library(rio)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
keytypes(org.Dr.eg.db)

data <- rio::import("all_gene.txt")

#GO富集分析
enrich.go <- enrichGO(gene = data$SYMBOL,  #待富集的基因列表
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

#write.table(enrich.go, 'all_gene_enrich.go.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#KEGG富集分析
gene <- bitr(data$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dr.eg.db")
enrich.kegg <- enrichKEGG(gene = gene$ENTREZID,,
                          organism = "dre",
                          keyType = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)


df_kegg <- as.data.frame(enrich.kegg)
#write.table(enrich.kegg, 'all_gene_enrich.kegg.txt', sep = '\t', row.names = FALSE, quote = FALSE)


library(clusterProfiler)
library(rio)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
keytypes(org.Dr.eg.db)

data <- rio::import("sig_gene.txt")

#GO富集分析
enrich.go <- enrichGO(gene = data$SYMBOL,  #待富集的基因列表
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'MF',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

ego_data_frame <- as.data.frame(enrich.go)

rich_fold <- function(x){
  generatio <- str_split(x[3],"/",simplify = T)
  bgratio <- str_split(x[4],"/",simplify = T)
  richratio <- (as.numeric(generatio[1])*as.numeric(bgratio[2]))/
    (as.numeric(generatio[2])*as.numeric(bgratio[1]))
  return(richratio)
}

ego_data_frame$Rich_fold <- apply(ego_data_frame,1,rich_fold)

# ego_data_frame <- ego_data_frame %>% 
#   top_n(10,Rich_fold)

dotplot(enrich.go,font.size=12,color="p.adjust",showCategory = 15,
        orderBy="Rich_fold",title="Enrichment GO MF Top15")

ggsave("sig_gene_enrich_go_MF.pdf",width=20,height=25,units = "cm")

#write.table(enrich.go, 'sig_gene_enrich.go.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#KEGG富集分析
gene <- bitr(data$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dr.eg.db")
enrich.kegg <- enrichKEGG(gene = gene$ENTREZID,,
                          organism = "dre",
                          keyType = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)


df_kegg <- as.data.frame(enrich.kegg)


#write.table(enrich.kegg, 'sig_gene_enrich.kegg.txt', sep = '\t', row.names = FALSE, quote = FALSE)

########可视化
library(readxl)
data <- read_excel("view.xlsx",sheet = 4)
data$Description <- factor(data$Description, levels = data$Description)

p <- ggplot(data=data, aes(x=Description, y=-log10(p.adjust), fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  theme_test() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 31)) +  # 移除纵轴的扩展空间
  scale_x_discrete(labels=data$Description) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched Pathway Terms")
p

ggsave(filename = "all_gene_go_padj_top10.pdf", width = 10, height = 6)


########可视化
library(readxl)
data <- read_excel("view.xlsx",sheet = 4)
data$Description <- factor(data$Description, levels = data$Description)

p <- ggplot(data=data, aes(x=Description, y=-log10(p.adjust), fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  theme_test() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 31)) +  # 移除纵轴的扩展空间
  scale_x_discrete(labels=data$Description) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched Pathway Terms")
p

ggsave(filename = "all_gene_go_padj_top10.pdf", width = 10, height = 6)


# 可视化KEGG
library(readxl)
data <- read_excel("view.xlsx",sheet = 5)

data <- data %>%
  mutate(GeneRatioValue = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) %>% 
  mutate(BgRatioValue = as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio))) %>%
  mutate(enrichRatio = GeneRatioValue / BgRatioValue) %>% 
  arrange(desc(enrichRatio)) %>% .[1:15,]

data$Description <- factor(data$Description, levels = data$Description)

p <- ggplot(data=data, aes(x=Description, y=enrichRatio, fill=p.adjust)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +
  theme_test() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +  # 移除纵轴的扩展空间
  scale_x_discrete(labels=data$Description) +
  xlab("KEGG term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched Pathway Terms")  +
  scale_fill_gradient(low = "#CAD3C8", high = "#2F4550")
p

ggsave(filename = "all_gene_kegg_padj_top10.pdf", width = 10, height = 6)


rich_fold <- function(x){
  generatio <- str_split(x[3],"/",simplify = T)
  bgratio <- str_split(x[4],"/",simplify = T)
  richratio <- (as.numeric(generatio[1])*as.numeric(bgratio[2]))/
    (as.numeric(generatio[2])*as.numeric(bgratio[1]))
  return(richratio)
}

library(stringr)
data$Rich_fold <- apply(data,1,rich_fold)

dotplot(data,font.size=8,color="p.adjust",showCategory = 15,
        x = sort(data$Rich_fold,decreasing=T))  # 画气泡图



library(ggpie)
library(ggplot2)

data <- rio::import("element.txt")
colnames(data) <- c("group","count")
"Intergenic" "CDS"        "3' UTR"     "Intron"     "5' UTR"
data$group <- factor(data$group, levels = c("3' UTR","Intron","5' UTR","Intergenic","CDS"))

# 使用ggpie绘图，并旋转图形30度
ggpie(data = data , group_key = "group", count_type = "count",
      label_info = "ratio", label_type = "circle",
      label_size = 4, label_pos = "out",border_size = 0.25,border_color = "white")
  

ggsave(filename = "distribution.pdf",height = 5,width = 5)



#20240505晚上

library(clusterProfiler)
library(rio)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
keytypes(org.Dr.eg.db)

data <- rio::import("sig_gene.txt")

#GO富集分析
enrich.go <- enrichGO(gene = data$SYMBOL,  #待富集的基因列表
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'MF',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

ego_data_frame <- as.data.frame(enrich.go)

rich_fold <- function(x){
  generatio <- str_split(x[3],"/",simplify = T)
  bgratio <- str_split(x[4],"/",simplify = T)
  richratio <- (as.numeric(generatio[1])*as.numeric(bgratio[2]))/
    (as.numeric(generatio[2])*as.numeric(bgratio[1]))
  return(richratio)
}

ego_data_frame$Rich_fold <- apply(ego_data_frame,1,rich_fold)

# ego_data_frame <- ego_data_frame %>% 
#   top_n(10,Rich_fold)

dotplot(enrich.go,font.size=12,color="p.adjust",showCategory = 15,
        orderBy="Rich_fold",title="Enrichment GO MF Top15")

ggsave("sig_gene_enrich_go_MF.pdf",width=20,height=25,units = "cm")



####################################################################################################


library(clusterProfiler)
library(rio)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
keytypes(org.Dr.eg.db)

data <- rio::import("sig_gene.txt")

#GO富集分析
enrich.go <- enrichGO(gene = data$SYMBOL,  #待富集的基因列表
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'CC',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

top_enrich_go <- enrich.go@result %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))

# 准备数据，将Description列的文本进行换行处理，限制每行的宽度
top_enrich_go$Description <- str_wrap(top_enrich_go$Description, width = 40)

# 绘制点图，并添加边框
p <- ggplot(top_enrich_go, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust), shape = 16) +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  # 使用蓝红色渐变
  scale_size(range = c(3, 10)) +  # 控制点的大小范围
  labs(title = "Enrichment GO CC Top15", x = "Gene Ratio", y = "GO Term") +
  theme_minimal() +  # 使用简洁的主题
  theme(text = element_text(size = 20),  # 设置文字大小
        legend.position = "right",  # 调整图例位置
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑色边框
        axis.text.y = element_text(size = 20))  # 设置y轴文字大小

p

ggsave(plot = p,"sig_gene_enrich_go_CC.pdf",width=25,height=25,units = "cm")




####################################################################################################


library(clusterProfiler)
library(rio)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
keytypes(org.Dr.eg.db)

data <- rio::import("sig_gene.txt")

#GO富集分析
enrich.go <- enrichGO(gene = data$SYMBOL,  #待富集的基因列表
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'BP',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

top_enrich_go <- enrich.go@result %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))

# 准备数据，将Description列的文本进行换行处理，限制每行的宽度
top_enrich_go$Description <- str_wrap(top_enrich_go$Description, width = 40)

# 绘制点图，并添加边框
p <- ggplot(top_enrich_go, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust), shape = 16) +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  # 使用蓝红色渐变
  scale_size(range = c(3, 10)) +  # 控制点的大小范围
  labs(title = "Enrichment GO BP Top15", x = "Gene Ratio", y = "GO Term") +
  theme_minimal() +  # 使用简洁的主题
  theme(text = element_text(size = 20),  # 设置文字大小
        legend.position = "right",  # 调整图例位置
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑色边框
        axis.text.y = element_text(size = 20))  # 设置y轴文字大小

p

ggsave(plot = p,"sig_gene_enrich_go_BP.pdf",width=25,height=25,units = "cm")


####################################################################################################

data <- rio::import("sig_gene.txt")

gene <- bitr(data$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Dr.eg.db")
enrich.kegg <- enrichKEGG(gene = gene$ENTREZID,,
                          organism = "dre",
                          keyType = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

top_enrich_kegg <- enrich.kegg@result %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))

# 准备数据，将Description列的文本进行换行处理，限制每行的宽度
top_enrich_kegg$Description <- str_wrap(top_enrich_kegg$Description, width = 40)

# 绘制点图，并添加边框
p <- ggplot(top_enrich_kegg, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust), shape = 16) +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  # 使用蓝红色渐变
  scale_size(range = c(3, 10)) +  # 控制点的大小范围
  labs(title = "Enrichment KEGG Top15", x = "Gene Ratio", y = "GO Term") +
  theme_minimal() +  # 使用简洁的主题
  theme(text = element_text(size = 20),  # 设置文字大小
        legend.position = "right",  # 调整图例位置
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑色边框
        axis.text.y = element_text(size = 20))  # 设置y轴文字大小

p

ggsave(plot = p,"sig_gene_enrich_KEGG.pdf",width=25,height=25,units = "cm")

####################################################################################################
library("rtracklayer")
gtf_data = import('C:/Users/11356/Desktop/GRCz11.gtf') #gtf的路径
#这里使用import导入gtf文件， 生成一个GRangs对象
gtf_data = as.data.frame(gtf_data)

gtf_anno <- gtf_data %>% filter(type %in% c("exon","five_prime_utr","three_prime_utr")) %>% 
  select(seqnames, start, end, gene_id, gene_name, type) %>% 
  as.data.frame()

write.table(gtf_anno, "gtf_anno.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



# 加载必要的库
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

# 读入GTF文件
gtf_path <- "C:/Users/11356/Desktop/GRCz11.gtf"  # 修改为你的GTF文件路径
txdb <- makeTxDbFromGFF(gtf_path, format="gtf")

# 提取并选择每个基因的最长转录本
transcripts <- transcriptsBy(txdb, by="gene")
longest_transcripts <- lapply(transcripts, function(x) x[which.max(width(x))])
longest_tx <- do.call("c", longest_transcripts)

# 提取CDS, 5'UTR和3'UTR
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
five_utrs <- fiveUTRsByTranscript(txdb, use.names=TRUE)
three_utrs <- threeUTRsByTranscript(txdb, use.names=TRUE)




cat gtf_anno.txt | awk -v OFS="\t" -v FS="\t" \
'{print "chr"$1,$2-1,$3,$4,$5,$6}' > gtf_anno.bed

#bed文件取交集,只输出两个文件的交集
bedtools intersect -a promoter.bed -b gtf_anno.bed -wa -wb > intersect.bed


bedtools intersect -a gtf_anno.bed -b promoter.bed > intersect.bed



####################################################################################################
####################################################################################################


library(clusterProfiler)
library(rio)
library(org.Dr.eg.db)
library(ggplot2)
library(dplyr)
keytypes(org.Dr.eg.db)

data <- rio::import("sig_gene.txt")

#GO富集分析
enrich.go <- enrichGO(gene = data$SYMBOL,  #待富集的基因列表
                      OrgDb = 'org.Dr.eg.db',
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'all',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)

top_enrich_go <- enrich.go@result %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))

# 准备数据，将Description列的文本进行换行处理，限制每行的宽度
top_enrich_go$Description <- str_wrap(top_enrich_go$Description, width = 70)

# 绘制点图，并添加边框
p <- ggplot(top_enrich_go, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust), shape = 16) +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  # 使用蓝红色渐变
  scale_size(range = c(3, 10)) +  # 控制点的大小范围
  labs(title = "Enrichment GO Top20", x = "Gene Ratio", y = "GO Term") +
  theme_minimal() +  # 使用简洁的主题
  theme(text = element_text(size = 25),  # 设置文字大小
        legend.position = "right",  # 调整图例位置
        panel.border = element_rect(colour = "black", fill = NA, size = 0.05),  # 增加边框粗细
        panel.grid.major = element_line(size = 0.05, colour = "gray"),  # 主要网格线粗细和颜色
        panel.grid.minor = element_line(size = 0.05, colour = "gray"),  # 次要网格线粗细和颜色
        axis.text.y = element_text(size = 12),  # 设置y轴文字大小
        plot.margin = unit(c(1, 1, 1, 1), "cm"))  # 调整图形边距

p

ggsave(plot = p,"sig_gene_enrich_go_top20_test15.pdf",width=25,height=25,units = "cm")


