library(ggplot2)
library(reshape2)

Time <- c(0, 2, 3, 4, 6, 8)
LOD <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
High <- c(0, 1.05, NA, 2.1, 2.7, 2.8)
Medium <- c(0, 0.7, NA, 1.4, 2.1, 2.45)
Low <- c(0, 0.35, NA, 0.7, 1.05, 1.40)
Noise <- c(NA, NA, 0, 0.35, 1.05, 1.75)

df <- data.frame(Time, LOD, High, Medium, Low, Noise)

df1 <- melt(df, na.rm = FALSE, value.name = "Value", id = "Time", variable = "Variable")

ggplot(df, aes(x = Time)) +
  geom_line(aes(y = LOD), linetype = "dotdash") +
  geom_smooth(aes(y = High), method = lm, formula = y ~ poly(x, 4), se = FALSE, color = "navy") +
  geom_smooth(aes(y = Medium), method = lm, formula = y ~ poly(x, 4), se = FALSE, color = "blue") +
  geom_smooth(aes(y = Low), method = lm, color = "steelblue") +
  geom_smooth(aes(y = Noise), method = lm, se = FALSE, color = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.01)) +
  theme_classic() +
  scale_colour_manual(name="legend", values=c("navy", "blue", "steelblue", "red")) +
  ylab("Metabolite relative abundance") +
  xlab("Time (hours)") +
  ggtitle("Exudates sampling time model")

ggsave(filename = "Exudates sampling time model.pdf", 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 20, height = 10)



ggplot(df1, aes(x = Time, y = Value)) +
  geom_smooth(aes(colour = Variable), method = lm, formula = y ~ poly(x, 3), se = FALSE) +
  scale_colour_manual(name="Legend", values=c("black", "navy", "blue", "steelblue", "red")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.01)) +
  theme_classic() +
  ylab("Metabolite relative abundance") +
  xlab("Time (hours)") +
  ggtitle("Exudates sampling time model")

ggsave(filename = "Exudates sampling time model-Legend.pdf", 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 20, height = 10)
