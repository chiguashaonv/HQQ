#=============== 关于多组之间数据的比较 ====================

#===== first step: load data =====
#这一步是为了把数据读取到文件里，read.csv就是读取文件后缀为csv的文件，所以你可以先把xlsx的文件在excel中保存成csv后缀的，csv的文件的意思是用逗号分割，所以如果你的数据中本身就有逗号，就不要保存成csv，而是保存成txt
#这一步中你自己做需要改动的也就是文件的地址，也就是[~/project/other.data/hqq/]这部分，把这部分改成你自己的地址就行了
gene <- read.csv("~/project/other.data/hqq/rnaseq.csv",header=T,stringsAsFactors = F)

#读取进来之后我们先看一下数据的样子，发现第一列是基因名，一般画图我们都会把基因名字设置为行名
rownames(gene) <- gene[,1] #rownames:设置gene这个数据框的行名，gene[,1]表示gene这个数据框的第一列，第二列就是gene[,2]
#然后删除掉gene的第一列
gene <- gene[,-1] #gene[,-1]就代表删除第一列，同样的，gene[,-2]就代表删除第二列

#这个时候，为了处理的方便，我们需要把数据的行列颠倒一下
gene.t <- as.data.frame(t(gene)) #t 这个命令就是颠倒行列的，但是数据框颠倒之后会变成matrix，所以我们用as.data.frame这个命令就是把转换之后的数据重新变成数据框
#数据框 data.frame和矩阵matrix的差别是什么？ 当我们对数据进行处理的时候，如果有一个矩阵，名字叫gene，列名包括“cell",”name",那么我们取第一列我们只能通过gene[,1]来达到目的
#但是，如果是data.frame，那么我们既可以通过gene[,1]也可以通过gene$cell来取出第一列

#接下来，我们需要对一个基因的每一个数据进行一个命名？就是说，这个数值，他是哪个细胞里的数值
# gene.t$celltype <- unlist(sapply(rownames(gene.t),function(x) substr(x,1,nchar(x)-1) ))#这个命令其实是R里面一个特殊的循环，现阶段可以试着理解,下面的这行命令是一样的目的
gene.t$celltype <- c(rep("active",3),rep("IL.15",3),rep("IL.2",3),rep("NA",3),rep("Q2",3),rep("S21",3))#rep这个命令呢，就是重复，rep("active",3)就是把active重复三遍

#===== 差异分析 =====
#anova
#举个具体基因的例子
gys1 <- aov(gene.t$Gys1~gene.t$celltype) #aov就是单因素方差分析，我们这个命令的目的就是看看celltype对于基因的表达有没有影响，有没有显著性
#aov的通常用法都是aov(formula),formula的格式就是y~x,y是因变量，基因表达放在Y,x是自变量，自变量可以有好多个，这里我们只研究细胞类型，
#比如说再加上浓度，那么命令可以这么写fit <- aov(gene.exp~celltype+conc),这个命令的目的就是研究这两个变量对于基因表达的影响

#查看细胞类型对于基因表达的影响是否显著
summary(gys1)
#==============================================================
#                Df Sum Sq Mean Sq F value   Pr(>F)
#gene.t$celltype  5 2687.7   537.5   27.95 3.23e-06 ***
#Residuals       12  230.8    19.2
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#这边我放一个结果解释一下，Pr(>F)这列，指的就是P值，我们可以看出，对于gys1这个基因来说，细胞类型的影响很大

#那么怎么看两组两组之间比较呢，很简单，一个命令就足够了
p.value <- TukeyHSD(gys1)
#还可以画出来这个结果
plot(TukeyHSD(gys1))#图上的线，是每两组的平均值的差值的置信区间，只要这条线没有经过0,那么这两组之间就是显著的，因为这个比较的组太多，所以纵坐标可能显示的不完全，这个纵坐标和上一个命令的结果的行名是完全对应的

