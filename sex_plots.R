
# NOTE set wd

# read the CSVs
x_depth <- read.csv("x_depth.csv", header = FALSE)
y_depth <- read.csv("y_depth.csv", header = FALSE)

# give the columns names to make it easier
colnames(x_depth) <- c("Sample", "X_Depth")
colnames(y_depth) <- c("Sample", "Y_Depth")

# sanity 
head(x_depth)
head(y_depth)

# merge on the sample column
depths <- merge(x_depth, y_depth, by = "Sample")

# sanity
head(depths)

# make the graphs

# load ggplot 
library(ggplot2)

# x depth plot
ggplot(depths, aes(x = Sample, y = X_Depth)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4)) +
  labs(title = "X Chromosome Depth per Sample", x = "Sample", y = "X Depth")

# y depth plot
ggplot(depths, aes(x = Sample, y = Y_Depth)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 4)) +
  labs(title = "Y Chromosome Depth per Sample", x = "Sample", y = "Y Depth")

