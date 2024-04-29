setwd("C:/Users/meeth/Desktop/RP2Files")
tbl <- read.csv('Canonicalgenomecoverage.csv')
tbl2 <- read.csv('Canonicalexomecoverage.csv')

coverage_values = tbl$CDS_coverage
coverage_values2 = tbl2$CDS_coverage



# Define bin edges and labels
bin_edges <- seq(1, 151, by = 10)
bin_labels <- paste(bin_edges[-length(bin_edges)], "-", bin_edges[-1])

# Create bins and assign coverage values to bins
coverage_bins <- cut(coverage_values, breaks = bin_edges, labels = bin_labels, right = FALSE)
coverage_bins2 <- cut(coverage_values2, breaks = bin_edges, labels = bin_labels, right = FALSE)

# Print summary of coverage bins
print(table(coverage_bins))
print(table(coverage_bins2))
hist(coverage_values,50)












