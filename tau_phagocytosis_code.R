library(ggplot2)

tau <- read.csv("/Users/gretchenclark/Documents/data analytics/tau_trip_no blanks.csv")
tau_3 <- tau
tau_3$image=NULL
tau_3$cell.number=NULL
tau_3$replicate=NULL
tau_3$mean=NULL
tau_3$min=NULL
tau_3 #get rid of extraneous columns
tau_4 <- na.omit(tau_3) #data without NA's

#boxplot
means <- aggregate(tau_4$background.subtracted.mean ~ tau_4$time.point, tau_4, mean) #aggregate means for boxplot
p <-ggplot(tau_4, aes(tau_4$time.point, tau_4$background.subtracted.mean, group = tau_4$time.point))

p + geom_jitter(alpha=.5) +
  geom_boxplot(alpha=.6, outlier.color = NA, color = "mediumpurple4") +
  xlab("\nTime Point (Hours)") + ylab("Mean Cell Fluorescence\n") + 
  theme_gray(base_size = 15) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3) #code for boxplot 

#confusion matrix
library(corrplot)
mcor1 <- cor(tau_4) #correlation plot of variables 
corrplot(mcor1)

#density plots 

d <- density(tau_4$background.subtracted.mean)
plot(d, main="density of mean")
polygon(d, col="pink", border="blue")

d2 <- density(tau_4$area)
plot(d2, main="density of area")
polygon(d2, col="orange", border="blue")


#PCA
library(factoextra)
pcatau <- prcomp(tau_4, scale = TRUE) #pca 
fviz_eig(pcatau)#scree plot
fviz_pca_ind(pcatau,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
) #scatterplot
fviz_pca_var(pcatau,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)#biplot


#PCA FROM http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

#ctree model
time <- as.factor(tau_4$time.point) #change timepoints to factor vector instead of numerical 
tau_5 <- data.frame(time, area, max, background.subtracted.mean) #make data frame to work off of 
tau_5 #check to make sure column names and observations are right 
library(ISLR) 
library(party) #ctree package
install.packages("partykit")
library(partykit) #ctree control package 
t4 <- ctree_control(mincriterion = 0.999, minsplit = 6, maxdepth = 3) #control command for tree 
tr5 <- ctree(time ~., data =tau_5, control = t4) #ctree modelk with pre-pruning
plot(tr5) #plotted ctree 

