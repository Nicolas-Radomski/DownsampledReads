#### Downsampling comparison ####

# clean environment
rm(list=ls())

# set working directory
setwd("/PathTo/RstudioWorkingDirectory/DownsampledReads")

# install regular packages
install.packages("ggplot2")
install.packages("plyr")
install.packages("reshape2")
install.packages("ggpubr")
install.packages("ggpmisc")

# call library
library(ggplot2)
library(plyr)
library(reshape2)
library(ggpubr)
library(ggpmisc)

# read dataframe
data_long = read.table("Listeria-cgMLST-Additional-file-1.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")

# check dimension
dim(data_long)
# => [1] 840   7

# check 10 first lines
head(data_long, 10)

# check nature of variables (integer or factor)
str(data_long)

# change interger as factor for variables targeted_kmer_depth and targeted_read_depth
data_long$targeted_kmer_depth = as.factor(data_long$targeted_kmer_depth)
data_long$targeted_read_depth = as.factor(data_long$targeted_read_depth)

# check nature of variables (integer or factor)
str(data_long)

# reorganize levels variables
levels(data_long$targeted_read_depth)
data_long$targeted_read_depth = factor(data_long$targeted_read_depth, levels=c("100", "90", "80", "70", "60", "50", "40", "30", "20", "10"))
levels(data_long$targeted_read_depth)
levels(data_long$targeted_kmer_depth)
data_long$targeted_kmer_depth = factor(data_long$targeted_kmer_depth, levels=c("75", "68", "60", "53", "45", "38", "31", "23", "16", "8"))
levels(data_long$targeted_kmer_depth)

# rename levels of a factor
levels(data_long$targeted_read_depth)
data_long$targeted_depth = revalue(data_long$targeted_read_depth, c("100"="Dr100-Dk75","90"="Dr90-Dk68","80"="Dr80-Dk60","70"="Dr70-Dk53","60"="Dr60-Dk45","50"="Dr50-Dk38","40"="Dr40-Dk31","30"="Dr30-Dk23","20"="Dr20-Dk16","10"="Dr10-Dk8"))
levels(data_long$targeted_depth)

# rename variables
colnames(data_long)
names(data_long)[names(data_long) == "sample"] <- "Sample"
names(data_long)[names(data_long) == "targeted_depth"] <- "Targeted_depth"
names(data_long)[names(data_long) == "read_depth"] <- "Estimated_read_depth"
names(data_long)[names(data_long) == "mapping"] <- "Mapping"
colnames(data_long)

# plot a 4-way figures with boxplots
p = ggplot(data = data_long, aes(x = Targeted_depth, y = Estimated_read_depth)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 8, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Estimated read depth (X) from BBMap or INNUca", limits = c(0, 110), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X) from BBNorm") +
  theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
        strip.text.x = element_text(size=16, face = "bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  facet_grid(Mapping ~ reference_genome)
p
plot(p)
ggsave("downsampling-depth.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("downsampling-depth.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# pass from long to short dataframe
data_short <- dcast(data_long, formula = Sample+reference_genome+Targeted_depth~Mapping, value.var = "Estimated_read_depth")

# check variables
str(data_short)

# check dimensions
dim(data_short)
# => [1] 420   5

# rename variables
colnames(data_short)
names(data_short)[names(data_short) == "reference_genome"] <- "Reference"
colnames(data_short)

#plot linear correlations and perform linear regressions
my.formula <- y ~ x
p = ggplot(data = data_short, aes(x = BBMap, y = INNUca, group=Reference, color=Reference)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = Reference, color = Reference, size = Reference)) +
  geom_smooth(aes(shape = Reference, color = Reference), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF', '#008000')) +
  scale_size_manual(values=c(3,3,3)) +
  theme(legend.position="bottom") +
  scale_x_reverse(name = "Estimated read depth (X) from BBMap", limits = c(105, 10), breaks = c(100,90,80,70,60,50,40,30,20,10)) +
  scale_y_continuous(name = "Estimated read depth (X) from INNUca", limits = c(10, 105), breaks = c(10,20,30,40,50,60,70,80,90,100)) +
  stat_poly_eq(formula = my.formula, 
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
    parse = TRUE,label.y = "top", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6)
p
plot(p)
ggsave("downsampling-lines-C.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("downsampling-lines-C.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# mean and SD for breadth
## all ATCC included
str(data_long$read_breadth)
data_long_BBMap=subset(data_long,Mapping == "BBMap")
str(data_long_BBMap)
mean(data_long_BBMap$read_breadth)
# => 99.34363
sd(data_long_BBMap$read_breadth)
# => 0.07713936

## ATCC19114
data_long_BBMap_ATCC19114=subset(data_long_BBMap,reference_genome == "ATCC19114")
mean(data_long_BBMap_ATCC19114$read_breadth)
# => [1] 99.30392
sd(data_long_BBMap_ATCC19114$read_breadth)
# => [1] 0.06665401

## ATCC19115
data_long_BBMap_ATCC19115=subset(data_long_BBMap,reference_genome == "ATCC19115")
mean(data_long_BBMap_ATCC19115$read_breadth)
# => [1] 99.35668
sd(data_long_BBMap_ATCC19115$read_breadth)
# => [1] 0.06852175

## ATCCBAA679
data_long_BBMap_ATCCBAA679=subset(data_long_BBMap,reference_genome == "ATCCBAA679")
mean(data_long_BBMap_ATCCBAA679$read_breadth)
# => [1] 99.37028
sd(data_long_BBMap_ATCCBAA679$read_breadth)
# => [1] 0.07952565

# mean and SD for depth
## all mapper included
ddply(data_long_BBMap, .(targeted_read_depth), summarize,  depth=mean(Estimated_read_depth), breadth=sd(Estimated_read_depth))
comment <- scan(what="character")
targeted_read_depth     depth   breadth
1                  100 101.49934 1.7405643
2                   90  92.06038 1.3490967
3                   80  81.20483 1.1945943
4                   70  71.68152 1.1431326
5                   60  60.73363 1.0384165
6                   50  51.13929 0.9248908
7                   40  41.59631 0.7338545
8                   30  30.87276 0.4860090
9                   20  21.56162 0.2892141
10                  10  10.93101 0.1416479

rm(comment)

## BBMap and INNUca mappers included
str(data_short)
ddply(data_short, .(Targeted_depth), summarize, BBMap_mean=mean(BBMap), BBMap_sd=sd(BBMap), INNUca_mean=mean(INNUca), INNUca_sd=sd(INNUca))
comment <- scan(what="character")
Targeted_depth BBMap_mean  BBMap_sd INNUca_mean INNUca_sd
1      Dr100-Dk75  101.49934 1.7405643    97.34810 2.1103765
2       Dr90-Dk68   92.06038 1.3490967    88.57167 2.2439826
3       Dr80-Dk60   81.20483 1.1945943    78.94738 1.9789158
4       Dr70-Dk53   71.68152 1.1431326    69.20286 2.0535588
5       Dr60-Dk45   60.73363 1.0384165    58.69571 1.9389262
6       Dr50-Dk38   51.13929 0.9248908    49.32714 1.6889011
7       Dr40-Dk31   41.59631 0.7338545    40.05071 1.7267852
8       Dr30-Dk23   30.87276 0.4860090    29.74381 1.2500809
9       Dr20-Dk16   21.56162 0.2892141    21.15881 0.9643327
10       Dr10-Dk8   10.93101 0.1416479    11.23571 0.1411130

rm(comment)

## BBMap and INNUca mappers for ATCC19114
data_short_ATCC19114=subset(data_short,Reference == "ATCC19114")
ddply(data_short_ATCC19114, .(Targeted_depth), summarize, BBMap_mean=mean(BBMap), BBMap_sd=sd(BBMap), INNUca_mean=mean(INNUca), INNUca_sd=sd(INNUca))
comment <- scan(what="character")
Targeted_depth BBMap_mean  BBMap_sd INNUca_mean INNUca_sd
1      Dr100-Dk75  101.56426 1.6420479    98.21429 1.4646066
2       Dr90-Dk68   91.91806 1.4746304    89.29786 1.7392004
3       Dr80-Dk60   80.94335 1.2799037    78.71071 1.8111724
4       Dr70-Dk53   71.39818 1.1057507    69.50214 2.0035514
5       Dr60-Dk45   60.45139 0.9362476    58.68714 2.0129406
6       Dr50-Dk38   50.85777 0.8100870    49.50286 1.1788288
7       Dr40-Dk31   41.34457 0.6178873    40.10643 0.9563601
8       Dr30-Dk23   30.72116 0.3842883    29.87214 1.5611691
9       Dr20-Dk16   21.48291 0.2463212    20.94714 0.8959089
10       Dr10-Dk8   10.89034 0.1265280    11.27143 0.1204388

rm(comment)

## BBMap and INNUca mappers for ATCC19115
data_short_ATCC19115=subset(data_short,Reference == "ATCC19115")
ddply(data_short_ATCC19115, .(Targeted_depth), summarize, BBMap_mean=mean(BBMap), BBMap_sd=sd(BBMap), INNUca_mean=mean(INNUca), INNUca_sd=sd(INNUca))
comment <- scan(what="character")
Targeted_depth BBMap_mean  BBMap_sd INNUca_mean INNUca_sd
1      Dr100-Dk75  101.91241 1.4271331    96.20643 1.7512220
2       Dr90-Dk68   92.28549 1.3594681    87.24143 2.3385527
3       Dr80-Dk60   81.27268 1.2097439    78.25357 1.7086774
4       Dr70-Dk53   71.65799 1.0520474    67.75643 1.5990700
5       Dr60-Dk45   60.68190 0.8846797    58.14286 1.4678893
6       Dr50-Dk38   51.08552 0.7836481    48.57357 1.5710095
7       Dr40-Dk31   41.54696 0.6126870    38.69643 1.3189775
8       Dr30-Dk23   30.82356 0.3762439    29.53286 1.1123186
9       Dr20-Dk16   21.53121 0.2238438    20.73786 0.8668704
10       Dr10-Dk8   10.91991 0.1202191    11.20714 0.1206666

rm(comment)

## BBMap and INNUca mappers for ATCCBAA679
data_short_ATCCBAA679=subset(data_short,Reference == "ATCCBAA679")
ddply(data_short_ATCCBAA679, .(Targeted_depth), summarize, BBMap_mean=mean(BBMap), BBMap_sd=sd(BBMap), INNUca_mean=mean(INNUca), INNUca_sd=sd(INNUca))
comment <- scan(what="character")
Targeted_depth BBMap_mean  BBMap_sd INNUca_mean INNUca_sd
1      Dr100-Dk75  101.02134 2.0946098    97.62357 2.5575397
2       Dr90-Dk68   91.97759 1.2781785    89.17571 2.1262416
3       Dr80-Dk60   81.39846 1.1314207    79.87786 2.1540130
4       Dr70-Dk53   71.98839 1.2662298    70.35000 1.7210864
5       Dr60-Dk45   61.06761 1.2380687    59.25714 2.2318962
6       Dr50-Dk38   51.47458 1.1036298    49.90500 2.0331928
7       Dr40-Dk31   41.89741 0.8790945    41.34929 1.7190850
8       Dr30-Dk23   31.07356 0.6204311    29.82643 1.0867425
9       Dr20-Dk16   21.67072 0.3640566    21.79143 0.8418373
10       Dr10-Dk8   10.98279 0.1673812    11.22857 0.1772811

rm(comment)

# correlation test

## Pearson correlation coefficient
cor(x=data_short_ATCC19114$BBMap, y=data_short_ATCC19114$INNUca, method="pearson")
# => [1] 0.9987964

## p-value of correlation coefficient (Significance levels)
### ATCC19114
cor.test(x=data_short_ATCC19114$BBMap, y=data_short_ATCC19114$INNUca, method=c("pearson"))
comment <- scan(what="character")
Pearsons product-moment correlation
data:  data_short_ATCC19114$BBMap and data_short_ATCC19114$INNUca
t = 239.22, df = 138, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9983180 0.9991388
sample estimates: cor 0.9987964 

rm(comment)

### ATCC19115
cor.test(x=data_short_ATCC19115$BBMap, y=data_short_ATCC19115$INNUca, method=c("pearson"))
comment <- scan(what="character")
Pearsons product-moment correlation
data:  data_short_ATCC19115$BBMap and data_short_ATCC19115$INNUca
t = 224.06, df = 138, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  0.9980832 0.9990185
sample estimates: cor 0.9986283

rm(comment)

### ATCCBAA679
cor.test(x=data_short_ATCCBAA679$BBMap, y=data_short_ATCCBAA679$INNUca, method=c("pearson"))
comment <- scan(what="character")
Pearsons product-moment correlation
data:  data_short_ATCCBAA679$BBMap and data_short_ATCCBAA679$INNUca
t = 208.95, df = 138, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.9977969 0.9988718
sample estimates: cor 0.9984234 

rm(comment)
