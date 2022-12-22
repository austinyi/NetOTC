library(tidyverse)
library(stringr)
library(latex2exp)

homedir <- # redacted for anonymity
options(stringsAsFactors=FALSE)

## Read Data.
accuracy_df <- read.csv(file.path(homedir, "alignment_accuracies.txt"))
edge_accuracy_df <- read.csv(file.path(homedir, "edge_alignment_accuracies.txt"))

## Wrangle data
acc_df <- accuracy_df %>%
  group_by(Algorithm) %>%
  summarise("mean_acc"=mean(Accuracy),
            "max_acc"=max(Accuracy),
            "min_acc"=min(Accuracy),
            "sd_acc"=sd(Accuracy))
acc_df$Algorithm[acc_df$Algorithm == "OTC"] <- "NetOTC"
acc_df$Algorithm[acc_df$Algorithm == "OT-SD"] <- "OT"
acc_df$Algorithm <- factor(acc_df$Algorithm, levels = c("NetOTC","OT","GW","FGW","COPT"))

edge_acc_df <- edge_accuracy_df %>%
  group_by(Algorithm) %>%
  summarise("mean_acc"=mean(Accuracy),
            "max_acc"=max(Accuracy),
            "min_acc"=min(Accuracy),
            "sd_acc"=sd(Accuracy))
edge_acc_df$Algorithm[edge_acc_df$Algorithm == "OTC"] <- "NetOTC"
edge_acc_df$Algorithm[edge_acc_df$Algorithm == "OT-SD"] <- "OT"
edge_acc_df$Algorithm <- factor(edge_acc_df$Algorithm, levels = c("NetOTC","OT","GW","FGW","COPT"))

## Make plots
plot_id <- format(Sys.time(), "%m-%d-%y_%H-%M-%S")
### Vertex accuracy
ggplot(acc_df, aes(x=Algorithm, y=mean_acc, fill=factor(Algorithm))) + 
  geom_bar(stat="identity",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_acc-sd_acc, ymax=mean_acc+sd_acc), width=.3,
                position=position_dodge(.9)) +
  theme_minimal() + 
  theme(legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.title.align=0.5,
        legend.box.background = element_rect(colour="black"),
        plot.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=8),
        axis.title = element_text(size=12, face="bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ggtitle("Vertex Alignment Accuracy") + 
  ylim(0, 0.8) +
  labs(fill = "Algorithm",
       color = "Algorithm",
       y = "Accuracy",
       x= "") +
  geom_hline(yintercept=0.25, lwd=0.6, colour="black", linetype="dashed")

ggsave(file.path(homedir, paste0("vertex_accuracy_", plot_id, ".png")), width=5, height=5)

### Edge accuracy
ggplot(edge_acc_df, aes(x=Algorithm, y=mean_acc, fill=factor(Algorithm))) + 
  geom_bar(stat="identity",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_acc-sd_acc, ymax=mean_acc+sd_acc), width=.3,
                position=position_dodge(.9)) +
  theme_minimal() + 
  theme(legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.title.align=0.5,
        legend.box.background = element_rect(colour="black"),
        plot.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=8),
        axis.title = element_text(size=12, face="bold"),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ggtitle("Edge Alignment Accuracy") + 
  ylim(0, 0.8) +
  labs(fill = "Algorithm",
       color = "Algorithm",
       y = "Accuracy",
       x= "") +
  geom_hline(yintercept=1/16, lwd=0.6, colour="black", linetype="dashed")

ggsave(file.path(homedir, paste0("edge_accuracy_", plot_id, ".png")), width=5, height=5)

