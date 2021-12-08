#this is the code of my project report

dat <- read.csv("data.csv",header = T)

tomor <- dat[dat$Type=="Tumor",]
NAT <- dat[dat$Type=="NAT",]
row.names(NAT) <- seq(1,99)

tomor.pair <- tomor[tomor$Participant %in% NAT$Participant,]
row.names(tomor.pair) <- seq(1,99)

tumor_mutate <- apply(tomor.pair[,c(34:46)],1,sum)
NAT_mutate <- apply(NAT[,c(34:46)],1,sum)

RNA <- read.csv("RNA.csv",header = T)
RNA <- as.data.frame(t(RNA))
colnames(RNA) <- RNA[1,]


gene <- read.csv("gene.csv",header = T)
name <- gene[-1,1]

RNA_filter <- RNA[RNA$id %in% name,]
RNA_filter <- RNA_filter[,-c(2:4)]
row.names(RNA_filter) <- RNA_filter[,1]

#癌和癌旁的基因表达量
tumor_gene <- RNA_filter[,tomor.pair$Sample.ID]
NAT_gene <- RNA_filter[,colnames(RNA_filter) %in% NAT$Sample.ID]

noname <- c("C3L.02646", "C3N.03072", "C3N.03662", "C3N.03886", "C3N.04155")
tumor_gene <- tumor_gene[,!(colnames(tumor_gene) %in% noname)]

#删除含有缺失值的基因
complete_gene <- intersect(row.names(NAT_gene)[complete.cases(NAT_gene)],
                           row.names(tumor_gene)[complete.cases(tumor_gene)])
tumor_gene <- tumor_gene[complete_gene,]
NAT_gene <- NAT_gene[complete_gene,]

#转为数值型
NAT_gene <- as.data.frame(apply(NAT_gene,2,as.numeric))
row.names(NAT_gene) <- row.names(tumor_gene)
tumor_gene <- as.data.frame(apply(tumor_gene,2,as.numeric))
row.names(tumor_gene) <- row.names(NAT_gene)

#取总表达量作为某一个体的癌相关基因表达量
NAT_total_exp <- apply(NAT_gene,2,sum)
tumor_total_exp <- apply(tumor_gene,2,sum)

#检验两者表达是否有显著差异
#检验正态性看是否能T检验
shapiro.test(NAT_total_exp)
shapiro.test(tumor_total_exp)


#符号秩检验
wilcox.test(NAT_total_exp,tumor_total_exp,paired = T,alternative = "greater")

library(ggstatsplot)
gene_exp_comp <- data.frame(type=rep(c("NAT","tumor"),c(94,94)),
                            total_exp=c(NAT_total_exp,tumor_total_exp))

ggwithinstats(data = gene_exp_comp,
              x=type,
              y=total_exp,
              type = "nonparametric",
              title = "paired-comparision of total expression between NAT and tumor",
              ylab = "total expression",
              results.subtitle = F)


#癌相关基因总体表达量和基因突变数目的关系
tumor_gene_99sample <- RNA_filter[complete_gene,tomor.pair[,1]]
tumor_gene_99sample <- as.data.frame(apply(tumor_gene_99sample,2,as.numeric))
row.names(tumor_gene_99sample) <- complete_gene
tumor_total_exp_99sample <- apply(tumor_gene_99sample,2,sum)

plot(tumor_mutate~tumor_total_exp_99sample,
     main="dot plot of mutation geen numbers~total expression",
     xlab="tumor-associated geens' total expression",
     ylab="mutation geen numbers")

ggscatterstats(correlation_of_mut_num_exp,
               y=mutant_number,
               x=all_gene,
               title="dot plot of mutation geen numbers~total expression",
               xlab="tumor-associated geens' total expression",
               ylab="mutation geen numbers")


signif_gene <- c("TP53","PTEN","CDKN2A","KMT2D","NFE2L2",
                 "ARID1A","CUL3","BRCA2","KEAP1","SUZ12",
                 "NF1","PIK3CA","NOTCH1")

signif_gene_99sample <- tumor_gene_99sample[signif_gene,]
signif_total_exp_99sample <- apply(signif_gene_99sample,2,sum)
cor.test(signif_total_exp_99sample,tumor_mutate)

correlation_of_mut_num_exp <- data.frame(mutant_number=tumor_mutate,
                                         all_gene=tumor_total_exp_99sample,
                                         signif_gene=signif_total_exp_99sample)

cor.test(signif_total_exp_99sample,tumor_mutate,method = "spearman")
ggcorrmat(correlation_of_mut_num_exp,
          cor.vars.names = list("mutation number","all","significant"),
          title = "correlation test and comparison",
          type = "np")


#回过头再做一遍只考虑显著突变基因的癌变和正常配对样本表达的关系
NAT_siginf_gene <- NAT_gene[signif_gene,]
NAT_siginf_total_exp <- apply(NAT_siginf_gene,2,sum)
tumor_siginf_gene <- tumor_gene[signif_gene,]
tumor_siginf_total_exp <- apply(tumor_siginf_gene,2,sum)
wilcox.test(NAT_siginf_total_exp,tumor_siginf_total_exp)

siginif_gene_exp_comp <- data.frame(type=rep(c("NAT","tumor"),c(94,94)),
                                    total_exp=c(NAT_siginf_total_exp,tumor_siginf_total_exp))
ggwithinstats(data = siginif_gene_exp_comp,
              x=type,
              y=total_exp,
              type = "nonparametric",
              title = "paired-comparision of significant gene expression between NAT and tumor",
              ylab = "significant gene expression",
              results.subtitle = F)

wilcox.test(NAT_siginf_total_exp,tumor_siginf_total_exp,alternative = "less")


#找相关因素
analyse <- tomor[,c(7,10,11,12,14,16,17,18)]
row.names(analyse) <- tomor[,1]
analyse <- cbind(analyse,apply(tomor[,34:46],1,sum))
colnames(analyse)[9] <- "mutation number"
analyse$`mutant level` <- ifelse(analyse$`mutation number`>=2,"high","low")

#查看突变数目分布情况
mycolor <- brewer.pal(12,"Set3")
median(analyse[,9])
gghistostats(analyse,
             x=`mutation number`,
             binwidth = 1,
             type = "np",
             bin.args = list(fill=mycolor[1],alpha=0.4,col="white"),
             normal.curve = T,
             normal.curve.args = list(size=1,col="red"),
             title = "distribution of mutant genes of 108 sample")

#chisq检验+ANOVA
ggbarstats(analyse,
           y=Smoking.History_modified,
           x=`mutant level`,
           type = "np",
           title="mutant level in different somking history",
           xlab = "smoking history")


ggbetweenstats(analyse,
           x=Smoking.History_modified,
           y=`mutation number`,
           title="mutant level in different somking history",
           xlab = "smoking history",
           pairwise.comparisons = T,
           results.subtitle = F)

kruskal.test(`mutation number`~Smoking.History_modified,data = analyse)

ggbarstats(analyse,
           y=Country.of.Origin,
           x=`mutant level`,
           title="mutant level in different country",
           xlab = "country",
           proportion.test = T)

ggbetweenstats(analyse,
               x=Country.of.Origin,
               y=`mutation number`,
               title="mutant level in different country",
               xlab = "country",
               pairwise.comparisons = T,
               results.subtitle = F)

kruskal.test(`mutation number`~Country.of.Origin,data = analyse)
#和年龄的关系
gghistostats(analyse,
             x=Age,
             type = "np",
             bin.args = list(fill=mycolor[10],alpha=0.4,col="white"),
             normal.curve = T,
             normal.curve.args = list(size=1,col="red"),
             title = "distribution of age of 108 sample")
analyse$`age group` <- ifelse(analyse$Age>=67,"old","young")

ggbarstats(analyse,
           y=`age group`,
           x=`mutant level`,
           title="mutant level in different age",
           xlab = "age group",
           package = "RColorBrewer",
           palette = "Accent")

ggbetweenstats(analyse,
               x=`mutant level`,
               y=Age,
               title="age distribution in different mutant level",
               xlab = "mutant level",
               package = "RColorBrewer",
               palette = "Accent")
kruskal.test(`mutation number`~Country.of.Origin,data = analyse)

#其他因素
ggbarstats(analyse,
           y=Gender,
           x=`mutant level`,
           title="mutant level in different gender",
           xlab = "gender",
           package = "RColorBrewer",
           palette = "Pastel1")

ggbarstats(analyse,
           y=Ethnicity_mod,
           x=`mutant level`,
           title="mutant level in different ethnicity",
           xlab = "ethnicity",
           package = "RColorBrewer",
           palette = "Pastel1")


ggbarstats(analyse,
           y=Ethnicity_mod,
           x=`mutant level`,
           title="mutant level in different ethnicity",
           xlab = "ethnicity",
           package = "RColorBrewer",
           palette = "Pastel1")

analyse$`daily cigarettes` <- ifelse(analyse$Cigarettes.per.Day=="Unknown",NA,as.numeric(analyse$Cigarettes.per.Day))
analyse$`daily cigarettes group` <- analyse$`daily cigarettes`
for(i in 1:nrow(analyse)){
  if(is.na(analyse$`daily cigarettes`[i]))analyse$`daily cigarettes group`[i] <- NA
  else if(analyse$`daily cigarettes`[i]<20)analyse$`daily cigarettes group`[i] <- "low"
  else if(analyse$`daily cigarettes`[i]==20)analyse$`daily cigarettes group`[i] <- "median"
  else if(analyse$`daily cigarettes`[i]>20)analyse$`daily cigarettes group`[i] <- "high"
}

ggbarstats(analyse,
           y=`daily cigarettes group`,
           x=`mutant level`,
           title="mutant level in different daily-cigarettes group",
           xlab = "daily-cigarettes group",
           package = "RColorBrewer",
           palette = "Set1")

shapiro.test(analyse$Pack.Years.Smoked)
bartlett.test(Pack.Years.Smoked~`mutant level`,analyse)

gghistostats(analyse,
             x=Pack.Years.Smoked,
             bin.args = list(fill=mycolor[12],alpha=0.4,col="white"),
             normal.curve = T,
             normal.curve.args = list(size=1,col="red"),
             title = "distribution of yearly-cigarette packs of 108 sample")
ggbetweenstats(analyse,
               x=`mutant level`,
               y=Pack.Years.Smoked,
               title="mutant level in different yearly-cigarette packs",
               xlab = "mutant level",
               package = "RColorBrewer",
               palette = "Set1",
               ylab = "yearly-cigarette packs")

ggbarstats(analyse,
           y=Secondhand.Smoke,
           x=`mutant level`,
           title="mutant level in different secondhand-smoke group",
           xlab = "secondhand-smoke group",
           package = "RColorBrewer",
           palette = "Set1")

analyse$gender <- ifelse(analyse$Gender=="female",0,1)


model <- glm(`mutation number`~Age+Pack.Years.Smoked+`daily cigarettes`+
               gender+Age*gender,data = analyse,family = poisson)
ggcoefstats(model)


X <- analyse[,c("mutation number","Age","Pack.Years.Smoked",
                "daily cigarettes","gender")]
X <- X[complete.cases(X),]
XX <- cor(X)
kappa(XX,exact = T)

#逐步回归
model.all <- glm(`mutation number`~Age+Pack.Years.Smoked+`daily cigarettes`+
                   gender+Age*gender+Age*Pack.Years.Smoked+
                   Age*`daily cigarettes`+gender*`daily cigarettes`+
                   gender*Pack.Years.Smoked,data = X,family = poisson)
model.new <- step(model.all,direction = "both")


ggscatterstats(analyse,
               x=`daily cigarettes`,
               y=`mutation number`)

ggscatterstats(analyse,
               x=Pack.Years.Smoked,
               y=`mutation number`)

ggscatterstats(analyse,
               x=Age,
               y=`mutation number`)
