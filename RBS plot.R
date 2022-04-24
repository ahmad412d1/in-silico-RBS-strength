getwd()

#import library
library(dplyr)
library(rstatix)
library(plotly)
library(ggpubr)

#read file
df <- read.csv("../dataset_rbs.csv", header=TRUE)

#A. visualize the RBS strength using plotly
fig1 <- df %>%
  plot_ly(
    x = ~Gene_name,
    y = ~pos.27.TIR.LOG,
    split = ~RBS_name,
    type = "bar",
    box = list(
      visible = T
    ),
    menaline = list(
      visible = T
    )
  )

#box plot on translation rate 
fig1 <- fig1 %>%
  layout(
    xaxis = list(
      title = "Genes with various RBS"
      ,categoryorder = "category ascending",
      categoryarray = ~Gene_name)
    ,
    yaxis = list(
      title = "Translation rate (Log10)",
      zeroline = F
    )
  )

fig1

#dG vs translation rate
fig2 <- df %>%
  plot_ly(
    x = ~pos.27.dG_total,
    y = ~pos.27.TIR.LOG,
    split = ~RBS_name,
    type = "scatter",
    box = list(
      visible = T
    ),
    menaline = list(
      visible = T
    )
  )

#box plot on translation rate 
fig2 <- fig2 %>%
  layout(
    xaxis = list(
      title = "dG total"
      ,categoryorder = "array",
      categoryarray = ~pos.27.TIR)
    ,
    yaxis = list(
      title = "Translation rate",
      zeroline = T
    )
  )

fig2

#B. Simple regression
#https://www.scribbr.com/statistics/linear-regression-in-r/
##normality check
hist(df$pos.27.dG_total)
hist(df$pos.27.TIR.LOG)

##Linearity
plot(pos.27.dG_total ~ pos.27.TIR.LOG, data = df)

##linear regression analysis
dG.TIR.lm <- lm(pos.27.dG_total ~ pos.27.TIR.LOG, data = df)
summary(dG.TIR.lm)

##visualization
dG.TR.graph <- ggplot(df, aes(x=pos.27.dG_total, y=pos.27.TIR.LOG)) +
  geom_point() +
  geom_smooth(method = "lm", col="blue") +
  stat_regline_equation(label.x = -5, label.y = 5) +
  theme_bw() +
  labs(title = "Relation between dG total (energy) vs translation rate",
       x = "dG Total",
       y = "Translation rate (Log10)")

dG.TR.graph

#C. clustering
#https://www.statology.org/k-means-clustering-in-r/
library(cluster)
library(factoextra)
library(gridExtra)
library(tidyverse)

#drop columns that contain NA
#usefull web https://www.datasciencemadesimple.com/drop-variables-columns-r-using-dplyr/
df_m <- read.csv("../out_RBS-CDS_all calculations_for clustering.csv", header=TRUE, row.names = "ID")
df1 <- df_m[,!sapply(df_m, function(x) mean(is.na(x)))>0.3]

#drop "non numeric" columns and set RBS name as index
#drop https://www.statology.org/r-drop-column/
df2 <- select(df1,-c(1,2))
df3 <-  subset(df2,select = -c(pos.8.TIR, pos.8.ORF, pos.8.dG_total, pos.8.dG_mRNA_rRNA, 
              pos.8.dG_spacing, pos.8.dG_standby, pos.8.dG_start, pos.27.ORF, 
              pos.27.dG_standby, pos.27.dG_start))


#make cluster plots
cluster <- df3 %>% slice(1:33)
cluster

cluster_plot <- scale(cluster)
cluster_plot

#Kmean clustering
kmean <- kmeans(cluster_plot, centers = 3, nstart = 33)
kmean

#plot results of final k-means model
fviz_cluster(kmean, data = cluster_plot)

#find means of each cluster
aggregate(df3, by = list(cluster=kmean$cluster), mean)

#add cluster assignment to original data
final_data <- cbind(df3, cluster = kmean$cluster)
head(final_data)

#export to .csv
write.csv(final_data, file="kmeancluster_RBS.csv")

#D. Stat test
#useful link https://www.datanovia.com/en/lessons/how-to-do-a-t-test-in-r-calculation-and-reporting/
#stat t-Test for LOG pos.27.TIR
stat.test.LOG <- df %>% 
  t_test(pos.27.TIR.LOG ~ RBS_name, var.equal = FALSE) %>%
  add_significance()
stat.test.LOG

#stat test for the absolute valvue pos.27.TIR
stat.test.a <- df %>%
  t_test(pos.27.TIR ~ RBS_name, var.equal = FALSE) %>%
  add_significance()
stat.test.a
