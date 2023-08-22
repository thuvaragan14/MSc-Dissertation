def Allele_Frequency(filename, output):
    with open(filename, "r") as f:
        # Open the output file
        with open(output, "w") as out_f:
            # Write the header line to the output file
            out_f.write("SNPID\tREF\tALT\tREF_FREQ\tALT_FREQ\n")
        
            # Iterate over each line in the vcf file
            for line in f:
                # Skip header lines
                if line.startswith("#"):
                    continue
            
                # Split the line into fields
                fields = line.strip().split("\t")
            
                # Get the reference allele and alternate alleles
                ref = fields[3]
                alt = fields[4]
            
                # Get the genotype information for all samples taking only the first three characters e.g. 0/1
                genotypes = fields[9:]
                genotypes_first_three = [genotype[:3] for genotype in genotypes]
            
                # Calculate the frequency of each allele
                ref_count = 0
                alt_count = 0
                for allele in genotypes_first_three:
                    if allele == "0/0":
                        ref_count += 2
                    elif allele == "0/1" or allele == "1/0":
                        ref_count += 1
                        alt_count += 1
                    elif allele == "1/1":
                        alt_count += 2
            
                # Write the result to the output file
                out_f.write("{}\t{}\t{}\t{}\t{}\n".format(
                    fields[0]+fields[1], ref, alt, ref_count, alt_count
                ))
