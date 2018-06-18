#========================== heatmap =================================

#===== first step: load data =====
#这一步是为了把数据读取到文件里，read.csv就是读取文件后缀为csv的文件，所以你可以先把xlsx的文件在excel中保存成csv后缀的，csv的文件的意思是用逗号分割，所以如果你的数据中本身就有逗号，就不要保存成csv，而是保存成txt
#这一步中你自己做需要改动的也就是文件的地址，也就是[~/project/other.data/hqq/]这部分，把这部分改成你自己的地址就行了
gene <- read.csv("~/project/other.data/hqq/rnaseq.csv",header=T,stringsAsFactors = F)

#读取进来之后我们先看一下数据的样子，发现第一列是基因名，一般画图我们都会把基因名字设置为行名
rownames(gene) <- gene[,1] #rownames:设置gene这个数据框的行名，gene[,1]表示gene这个数据框的第一列，第二列就是gene[,2]
#然后删除掉gene的第一列
gene <- gene[,-1] #gene[,-1]就代表删除第一列，同样的，gene[,-2]就代表删除第二列

#===== second step: plot heatmap =====
#接下来就可以开始画图了
#画热图的时候，我们经常使用的一个包叫做pheatmap，这个包画出来的热图好看，而且可调节性比较大，比较灵活
#如果你已经安装了这个包，就直接加载
library(pheatmap)
#如果你没有安装，就先安装
install.packages("pheatmap")
#然后再加载

#接下来就是画图
annotation_col = data.frame("Cell.type"=factor(c(rep("active",3),rep("IL-15",3),rep("IL-2",3),rep("naive",3),rep("Q2",3),rep("S21",3)),labels=c("active","IL-15","IL-2","naive","Q2","S21")))
rownames(annotation_col) <- colnames(gene)
annotation_colors = list(Cell.type=c("active"="#A7414A","IL-15"="#282726",
	"IL-2"="#6A8A82","naive"="#A37C27","Q2"="#563838","S21"="#6465A5"))[1]
p.gene <- pheatmap(gene,#gene就是我们上述的数据框，这个数据框行是基因名，列是样本，所以你画出来的图也是行是基因，列是样本
	# color=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),#设置热图的颜色，从高到低分别是红色，黄色和蓝色,这是其中一种方法，
	color=colorRampPalette(c("blue","white","red"))(100),#这是另外一种设置颜色的方法， 这种方法灵活很多，你可以自己修改颜色，所以我们这里使用这种方法
	border_color=NA,#设置每个小格子的边框的颜色，这里我们设置为NA，也就是没有边框
	scale="row",#对数据以行为单位进行正态化，也可以用scale="col",以列为单位进行正态化，或者scale="none"，不进行正态化，一般来说，都是以一个基因为单位进行正态化，因为我们的行是基因，所以我们选择以行为单位
	cluster_cols=F,#对列进行聚类，这里我们选择否，也可以选择是cluster_cols=T
	cluster_rows=T,#与上面一样，这不过是针对行
	clustering_method="complete",#聚类的方法，这个是默认的，你还可以选择"ward.D", "ward.D2", "single", "complete", "average" , "mcquitty", "median" or "centroid"
	show_colnames=T,show_rownames=T,#显示行名和列名，我们选择显示
	fontsize=9, fontsize_row=10,fontsize_col=10,#设置字体大小，fontsize_row设置行名字体大小，fontsize_col设置列名
    annotation_col=annotation_col,
    annotation_colors = annotation_colors)