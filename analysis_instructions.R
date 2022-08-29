#Commented for instructions to first-time users

#install.packages("agricolae")
#install.packages("dplyr")
#install.packages("ggplot2")

#load the package
library(agricolae)
library(dplyr)
library(ggplot2)

#set working domain. telling R which folder to look into.
#replace it with your own by selecting the folder in the
#bottom right menu
setwd("C:/Users/44799/Desktop/xandrine!/babeR")

#create a datatable called data1 from a csv file w headers 
data1 <- read.csv(file = "DA14.csv", header = TRUE)

reunion <- subset(data1, Group == "Reunion")
isolated <- subset(data1, Group == "Isolated")
control <- subset(data1, Group == "Control")

reunion_avg <- mean(reunion$percent)
isolated_avg <- mean(isolated$percent)
control_avg <- mean(control$percent)



###############Tukey Test cFos_TH ~ Group

#saves a linear model analysis of two variables in data1
model1 <- lm(data1$cFos_TH ~ data1$Group)

#returns a summary of the linear model analysis
summary(model1)

#runs anova analysis on the lm analysis
anova(model1)

#some pretty Tukey test specific stuff here, not super relevant
#unless youre doing another one

a1 <- aov(data1$cFos_TH ~ data1$Group)
posthoc1 <- TukeyHSD(x=a1, 'data1$Group', conf.level=0.95)

HSD1 <- HSD.test(model1, 'data1$Group')

#########Tukey Test %activeTH ~ Group
model2 <- lm(data1$percent ~ data1$Group)
summary(model2)
anova(model2)

a2 <- aov(data1$percent ~ data1$Group)
posthoc2 <- TukeyHSD(x=a2, 'data1$Group', conf.level=0.95)

HSD2 <- HSD.test(model2, 'data1$Group')

##############bad boi histogram timeee

#brings up a guide on how to use the function 'boxplot'
#?boxplot

#makes a basic boxplot for two variables in data1
boxplot(data1$cFos_TH ~ data1$Group)

#boxplot, two variables, x label, y label, y range, colours
boxplot(data1$percent ~ data1$Group, xlab = "Group",
        ylab = "Active Dopamine Neurones (%)", ylim = c(0, 100),
        col = c("lightblue", "magenta", "turquoise"))


df <- data.frame(group=c("Reunion", "Isolated", "Control"),
                 percent=c(reunion_avg, isolated_avg, control_avg))
barplot1 <- 
  ggplot(data = df, aes(x = group, y = percent)) + geom_bar(stat="identity", width = 0.5) +
         xlab("Behavioural Group") + ylab("Percent Double Labelled TH & cFOS Neurones (%)") +
         geom_errorbar(aes(ymin= percent, ymax= c(47.34, 52.52, 14.56)), width= .3)


# create a dataset
time <- c("0-30min", "0-30min", "30-60min", "30-60min")
Group <- c("Isolation", "Reunion", "Isolation", "Reunion")
contact <- c(58.875, 71.25, 91.5, 120)
data <- data.frame(time,Group,contact)

# Grouped
ggplot(data, aes(fill=Group, y=contact, x=time)) + 
  geom_bar(position= position_dodge(), stat="identity") + 
  xlab("Time (mins)") + ylab("N of Contact") +
  geom_errorbar(aes(x = time, ymin= c(58.875, 71.25, 91.5, 120), 
  ymax= c(89.15, 92.90, 131.43, 120)), width= .25,
  position = position_dodge(.9)) +
  scale_fill_manual(values=c("gray80","gray40"))

  










