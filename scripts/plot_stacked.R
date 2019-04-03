library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(plotly)
#read in required packages
require(data.table)
library(Rmisc)



dat1 = read.csv("se_km_3_len.tsv", header = FALSE, sep=' ', colClasses=c("factor","numeric", "numeric","numeric","numeric", "numeric","numeric","numeric"))
#dat2 = read.csv("new_counts.tsv", header = TRUE, sep=' ', colClasses=c("factor", "numeric", "numeric","numeric","numeric"))
names(dat1)= c("ID","min","mid","max","TP","FP","FN","TN")
#dat$ID <- reorder(dat$ID, rowSums(dat[,2:5]))
data_subset <- subset(dat1, select = c("ID","TP","FP","FN","TN"))

d <- arrange(data_subset, TP+FN)
d$ID <- reorder(d$ID, rowSums(d[,c(2,4)]))

#d[,2:5] <- d[,2:5] * (matrix((nrow(d)-1):0) * (100/nrow(d)))
d$FN <- -1*d$FN
d$TN <- -1*d$TN

md <- melt(d, id=(c("ID")))
ggplot(data=md, aes(x=ID, y=value, fill=factor(variable, levels=c("FP","TP","TN", "FN")))  ) + 
	geom_bar(stat="identity") +
	scale_fill_manual(values = c("#00BFC4","#F8766D","#00BFC4", "#F8766D"), name="Type") +
	#scale_fill_manual(values = c("#CC79A7","#0072B2", "#E69F00","#56B4E9"), name="Type") +
	geom_hline(yintercept = 0) +
	theme(legend.position = "none")
	
ggplot(dat1, aes(x=FP, y=TP)) + geom_point()

