# Read segmental duplication data
segmental_duplications = []
with open("segmental_duplications_data.txt", "r") as segment_file:
    next(segment_file)  # skip header
    for line in segment_file:
        fields = line.strip().split()
        segmental_duplications.append((fields[1], int(fields[2]), int(fields[3])))

# Read genome and exome data
genome_genes = set()
exome_genes = set()
with open("genome_exome_data.txt", "r") as data_file:
    next(data_file)  # skip header
    for line in data_file:
        fields = line.strip().split()
        gene_name = fields[0]
        start = int(fields[-4])
        end = int(fields[-3])

        # Check if gene overlaps with segmental duplication
        for chrom, dup_start, dup_end in segmental_duplications:
            if fields[1] == chrom and (start <= dup_end and end >= dup_start):
                # Check if gene is present in both genome and exome
                if gene_name in genome_genes:
                    exome_genes.add(gene_name)
                else:
                    genome_genes.add(gene_name)

# Read disease-causing genes
disease_genes = set()
with open("disease_genes_data.txt", "r") as disease_file:
    next(disease_file)  # skip header
    for line in disease_file:
        gene_name = line.strip().split()[0]
        disease_genes.add(gene_name)

# Find genes present in both exome and genome, and also in segmental duplications, and are disease-causing
common_genes = exome_genes.intersection(genome_genes).intersection(disease_genes)

# Print the list of common genes
print("Genes that are disease-causing, present in both genome and exome, and fall in segmental duplication regions:")
for gene in common_genes:
    print(gene)
