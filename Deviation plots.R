setwd("C:/Users/meeth/Desktop/RP2()")
library(ggplot2)
read_csv = read.csv("C:/Users/meeth/Desktop/RP2Files/Canonicalexomecoverage.csv")
y=read_csv$CDS_coverage

#Normal Graph
data <- data.frame(index = 1:length(y), coverage = y)
threshold <- 20
deviation_index <- min(which(data$coverage >= threshold))
plot <- ggplot(data, aes(x = index, y = coverage)) +
  geom_line() +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_vline(xintercept = deviation_index, linetype = "dashed", color = "blue") +
  labs(title = "Genome Coverage Deviation Visualization",
       x = "List of Canonical Transcripts of Genome Coverage",
       y = "Exon Coverage Values",
       caption = "Red dashed line: Threshold\nBlue dashed line: Deviation point") +
  theme_minimal()
print(plot)


#Histogram
data1 <- data.frame(coverage = y)
# Create the histogram plot
hist_plot <- ggplot(data1, aes(x = coverage)) +
  geom_histogram(binwidth = 7, fill = "blue", color = "black") +  # Adjust binwidth as needed
  labs(title = "CDS Coverage Frequency Distribution of Exome coverage",
       x = "Coverage Value",
       y = "Frequency") +
  theme_minimal()

# Display the histogram plot
print(hist_plot)
