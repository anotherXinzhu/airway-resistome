library(ggplot2)
library(dplyr)
library(ggpmisc)

load("../_data/TotalARG.abund.shannon.RData")
head(totalARGIndices)

ggplot(totalARGIndices) +
  geom_point(aes(x=ARGabundance, y=ARGshannon), shape = 21, fill = "#488DC3",color="black", alpha=0.6) +
  geom_smooth(aes(x=ARGabundance, y=ARGshannon)) 

# Using loess fitting, it was found that after a certain threshold, an increase in ARG abundance is accompanied by a significant decrease in ARG diversity.
# Therefore, plot the relationship between abundance and diversity in two segments:
# Use the segmented package to find the threshold of x.

library(segmented)
lm_model <- lm(ARGshannon ~ ARGabundance, data = totalARGIndices) # First, fit a linear regression model.
seg_model <- segmented(lm_model, seg.Z = ~ARGabundance) # Use the segmented package to find the optimal breakpoint.
summary(seg_model) # Check the fitting results: the breakpoint is 9441.933.

x_cutoff = 9441.933
totalARGIndices$Group <- ifelse(totalARGIndices$ARGabundance > x_cutoff, "High ARG","Low ARG")

totalARGIndices %>%
  group_by(Group) %>%
  summarise(p = cor.test(ARGabundance, ARGshannon)$p.value,
            r = cor.test(ARGabundance, ARGshannon)$estimate) 


p0 <- ggplot(totalARGIndices) +
  geom_point(aes(x=ARGabundance, y=ARGshannon), shape = 21, fill = "#488DC3",color="black", alpha=0.6) +
  geom_smooth(aes(x=ARGabundance, y=ARGshannon, group = Group), method = "lm") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Total ARG abundance (Coverage/Gbp)") +
  ylab("Shannon index of ARG") +
  geom_vline(xintercept = x_cutoff, linetype="dashed", color="gray")
p0



p.up <- ggplot(totalARGIndices) +
  geom_density(aes(x=ARGabundance), fill="skyblue3", 
               color="black", alpha=0.6, position = 'identity') +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank())


p.right <-ggplot(totalARGIndices) +
  geom_density(aes(x=ARGshannon), fill="skyblue3", 
    color="black", alpha=0.6, position = 'identity') +
  theme_classic() +
  coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())



library(ggpubr)

ggarrange(ggarrange(p.up, ggplot(), widths = c(0.85,0.15)),
          ggarrange(p0, p.right, widths = c(0.85,0.15)),
          nrow = 2, heights = c(0.15, 0.85)) 
