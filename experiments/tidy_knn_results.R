library(tidyverse)

datadir <- #redacted for anonymity
experiment <- "Yeast\\738293_390896"
k <- 5

data = data.frame(Algorithm=c(), Accuracy=c())
for (f in list.files(file.path(datadir, experiment), pattern="*.txt")){
  temp <- read.csv(file.path(datadir, experiment, f), header=TRUE, row.names=1)
  algorithms = row.names(temp)
  temp <- temp[,k]
  data <- rbind(data, data.frame(Algorithm=algorithms, Accuracy=temp))
}

data %>%
  group_by(Algorithm) %>%
  summarise("mean_acc"=mean(Accuracy),
            "max_acc"=max(Accuracy),
            "min_acc"=min(Accuracy),
            "sd"=sd(Accuracy))