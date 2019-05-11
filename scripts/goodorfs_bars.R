library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(plotly)
#read in required packages
require(data.table)
library(Rmisc)

#create a list of the files from your target directory
file_directory = "data/predicts/se_km_3_all/"

file_list <- list.files(path=file_directory)

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- data.frame()

#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  temp_data <- fread(paste(file_directory,file_list[i], sep=''), stringsAsFactors = F)
  dataset <- rbindlist(list(dataset, temp_data), use.names = T) #for each iteration, bind the new data to the building dataset
}
names(dataset)= c("ID","min","mid","max","TP","FP","FN","TN")

#dat1 = read.csv("old_counts.tsv", header = TRUE, sep=' ', colClasses=c("factor", "numeric", "numeric","numeric","numeric"))
#dat2 = read.csv("new_counts.tsv", header = TRUE, sep=' ', colClasses=c("factor", "numeric", "numeric","numeric","numeric"))
#names(dat)= c("ID","TP","FP","FN","TN")
#dat$ID <- reorder(dat$ID, rowSums(dat[,2:5]))

data_subset <- subset(dataset, select = c("ID","TP","FP","FN","TN"))

d <- arrange(data_subset, TP+FN)
d$ID <- reorder(d$ID, rowSums(d[,c(2,4)]))

#d[,2:5] <- d[,2:5] * (matrix((nrow(d)-1):0) * (100/nrow(d)))
d$FN <- -1*d$FN
d$TN <- -1*d$TN
#md <- melt(d, id=(c("ID")))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 600))) # 3 rows, 5 columns

s <- d[1:100,]
md <- melt(s, id=(c("ID")))
p1 <- ggplot(data=md, aes(x=ID, y=value, fill=factor(variable, levels=c("FP","TP","TN", "FN")))  ) + 
	geom_bar(stat="identity") +
	scale_fill_manual(values = c("#00BFC4","#F8766D","#00BFC4", "#F8766D"), name="Type") +
	#scale_fill_manual(values = c("#CC79A7","#0072B2", "#E69F00","#56B4E9"), name="Type") +
	geom_hline(yintercept = 0) +
	theme(legend.position = "none") +
	#ylim(-max(s[,2]+s[,3],-s[,4]-s[,5]) , max(s[,2]+s[,3],-s[,4]-s[,5]) ) +
	coord_cartesian(ylim=c(-max(s[,2]+s[,3]), max(s[,2]+s[,3]))) +
	ylab("orfs")

s <- d[101:200,]
md <- melt(s, id=(c("ID")))
p2 <- ggplot(data=md, aes(x=ID, y=value, fill=factor(variable, levels=c("FP","TP","TN", "FN")))  ) + 
	geom_bar(stat="identity") +
	scale_fill_manual(values = c("#00BFC4","#F8766D","#00BFC4", "#F8766D"), name="Type") +
	#scale_fill_manual(values = c("#CC79A7","#0072B2", "#E69F00","#56B4E9"), name="Type") +
	geom_hline(yintercept = 0) +
	theme(legend.position = "none") +
	#ylim(-max(s[,2]+s[,3],-s[,4]-s[,5]) , max(s[,2]+s[,3],-s[,4]-s[,5]) ) +
	coord_cartesian(ylim=c(-max(s[,2]+s[,3]), max(s[,2]+s[,3]))) +
	theme(axis.title.y=element_blank()) #, axis.text.y=element_blank())
	
s <- d[201:600,]	
md <- melt(s, id=(c("ID")))
p3 <- ggplot(data=md, aes(x=ID, y=value, fill=factor(variable, levels=c("FP","TP","TN", "FN")))  ) + 
	geom_bar(stat="identity") +
	scale_fill_manual(values = c("#00BFC4","#F8766D","#00BFC4", "#F8766D"), name="Type") +
	#scale_fill_manual(values = c("#CC79A7","#0072B2", "#E69F00","#56B4E9"), name="Type") +
	geom_hline(yintercept = 0) +
	theme(legend.position = "none") +
	#ylim(-max(s[,2]+s[,3],-s[,4]-s[,5]) , max(s[,2]+s[,3],-s[,4]-s[,5]) )+
	coord_cartesian(ylim=c(-max(s[,2]+s[,3]), max(s[,2]+s[,3]))) +
	theme(axis.title.y=element_blank()) #, axis.text.y=element_blank())

print(p1, vp = vplayout(1, 1:100))
print(p2, vp = vplayout(1, 101:200))
print(p3, vp = vplayout(1, 201:600))

## single plot of all
data_subset <- subset(dataset, select = c("ID","TP","FP","FN","TN"))
d <- arrange(data_subset, TP+FN)
d$ID <- reorder(d$ID, rowSums(d[,c(2,4)]))
#d <- d[1:1000,]
#d[,2:5] <- d[,2:5] * (matrix((nrow(d)-1):0) * (100/nrow(d)))
d[,2:5] <- log10(d[,2:5]+1)
d$FN <- -1*d$FN
d$TN <- -1*d$TN
md <- melt(d, id=(c("ID")))
ggplot(data=md, aes(x=ID, y=value, fill=factor(variable, levels=c("FP","TP","TN", "FN")))  ) + 
	geom_bar(stat="identity") +
	scale_fill_manual(values = c("#00BFC4","#F8766D","#00BFC4", "#F8766D"), name="Type") +
	#scale_fill_manual(values = c("#CC79A7","#0072B2", "#E69F00","#56B4E9"), name="Type") +
	geom_hline(yintercept = 0) +
	theme(legend.position = "none") +
	scale_y_continuous(breaks = c(-3,-2,-1,-0.3,0.3,1,2,3), labels = c("1000","100","10","1","1","10", "100","1000"))




