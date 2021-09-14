
genome_length = 4411532 #length of the Mycobacterium tuberculosis genome

fwd_sam_file = input("Enter the path and filename for the .sam file where reads were mapped to the + strand\n")
rev_sam_file = input("Enter the path and filename for the .sam file where reads were mapped to the - strand\n")
output_gff_file = input("Enter the path and filename for the output .gff file\n")

f = open(fwd_sam_file,'r') #SAM file from fwd genome alignment

genome_coverage_minus = [0 for i in range(genome_length + 1)] #1 bigger than genome_length to account for python counting

for x in f:
    SAM_line = x.split('\t')
    if SAM_line[1] == '16':
        genome_coverage_minus[int(SAM_line[3])] += 1
f.close()



f = open(rev_sam_file,'r') #SAM file from rev genome alignment

genome_coverage_plus = [0 for i in range(genome_length + 1)]

for x in f:
    SAM_line = x.split('\t')
    if SAM_line[1] == '16':
        genome_coverage_plus[genome_length - int(SAM_line[3]) +1] += 1
f.close()




f = open(output_gff_file,'w') #output gff

for x in range(len(genome_coverage_plus)):
    if genome_coverage_plus[x] > 0:
        f.write('NA\tNA1\tNA2\t')
        f.write(str(x) + '\t')
        f.write(str(x) + '\t')
        f.write(str(genome_coverage_plus[x]) + '\t.\t.\t.\n')

for x in range(len(genome_coverage_minus)):
    if genome_coverage_minus[x] > 0:
        f.write('NA\tNA1\tNA2\t')
        f.write(str(x) + '\t')
        f.write(str(x) + '\t')
        f.write(str(-genome_coverage_minus[x]) + '\t.\t.\t.\n')
f.close()
