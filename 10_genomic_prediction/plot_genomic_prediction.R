###This script is used to plot the results of genomic prediction accuracy###

#1. Load the necessary packages

library("data.table")
library("ggplot2")
library("gridExtra")
library("RColorBrewer")
library("ggpubr")

#2. Set the plot color

colors <- brewer.pal(8, "Dark2")[c(1)]

#3. Set the working directory that contains the genomic prediction results

setwd('/home/userpc377/TFM_data/prediction/ten_repeats/')

#4. Load the results for all the traits and create the plots

gwi <- fread("grain_width_prediction_10fold.txt", header = TRUE)

p1 <- ggplot() +
  geom_boxplot(data = gwi, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

gle <- fread("grain_length_prediction_10fold.txt", header = TRUE)

p2 <- ggplot() +
  geom_boxplot(data = gle, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

gxp <- fread("grainsxpanicle_prediction_10fold.txt", header = TRUE)

p3 <- ggplot() +
  geom_boxplot(data = gxp, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

ple <- fread("panicle_length_prediction_10fold.txt", header = TRUE)

p4 <- ggplot() +
  geom_boxplot(data = ple, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

npa <- fread("n_panicles_prediction_10fold.txt", header = TRUE)

p5 <- ggplot() +
  geom_boxplot(data = npa, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

hei <- fread("height_prediction_10fold.txt", header = TRUE)

p6 <- ggplot() +
  geom_boxplot(data = hei, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1,color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

d2h <- fread("days2heading_prediction_10fold.txt", header = TRUE)

p7 <- ggplot() +
  geom_boxplot(data = d2h, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

ale <- fread("allelopathy_prediction_10fold.txt", header = TRUE)

p8 <- ggplot() +
  geom_boxplot(data = ale, aes(x = Type, y = Correlation), fill = colors[1], alpha = 1, color = "black", linewidth = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::comma) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

#5. Present the plots

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4)

#6. Clean the environment

rm(list = ls())
gc()
