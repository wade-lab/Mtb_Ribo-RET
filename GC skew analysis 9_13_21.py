
import datetime
import os
import sys
import scipy
from scipy import stats

date = str(datetime.datetime.now().month) + "_" + str(datetime.datetime.now().day) + "_" + str(datetime.datetime.now().year)
time = str(datetime.datetime.now().hour) + "h_" + str(datetime.datetime.now().minute) + "m"




##############################################################################################
#                                          FUNCTIONS                                         #
##############################################################################################



def revcomp(sequence):
    revseq = sequence[::-1]

    revseq = revseq.replace('A','k')
    revseq = revseq.replace('C','l')
    revseq = revseq.replace('G','m')
    revseq = revseq.replace('T','n')
    
    revseq = revseq.replace('k','T')
    revseq = revseq.replace('l','G')
    revseq = revseq.replace('m','C')
    revseq = revseq.replace('n','A')

    return revseq

##############################################################################################


filepath = input("Enter the filepath where you want the output directory to be written\n")
filepath_input_files = input("Enter the filepath where the input files are located\n")
number_of_tails = int(input("One- or two-tailed Fisher's exact test? ['1' or '2']\n"))
if number_of_tails != 1 and number_of_tails != 2:
    print("Number must be 1 or 2. Please re-run")
    sys.exit()

filepath_output_files = filepath + 'Output_files_Overlapping_ORF_GC_skew_analysis_' + date + '_' + time + '/'
os.makedirs(filepath_output_files)

genome_fasta = input("Enter the name of the genome fasta file\n")
f = open(filepath_input_files + genome_fasta, 'r')
qqq = f.readline() #skip the header line
gen = 'x' + f.readline().split()[0].upper() #add an "x" at the beginning so the +1 genome position is position [1]
f.close()
revcompgen = 'x' + revcomp(gen) #revcomp genome sequence


gene_gff = input("Enter the name of the gff file with annotated genes to determine overlap\n")


ORF_file = input("Enter the name of the file with the list of ORFs to be analyzed\n") #strand, start codon coordinate, stop codon coordinate

output_filename = input("Enter name for output file ('.txt' will be added automatically)\n")



len_of_genome = 4411532 #length of the Mtb H37Rv genome





genome_coverage = [0 for x in range(len_of_genome)] #set all genome coordinates to 0

with open(filepath_input_files + gene_gff, 'r') as annotated_genes_gff: #any genome coordinate overlapping a gene is set to 1
    for x in annotated_genes_gff:
        for y in range(int(x.split()[3]), int(x.split()[4])):
            genome_coverage[y] = 1


ORF_list = [] #list of the ORFs, read in from file

with open(filepath_input_files + ORF_file, 'r') as ORFs:
    for x in ORFs:
        ORF_list.append([x.split()[0], int(x.split()[1]), int(x.split()[2])]) #list of lists: strand, start, stop
            

for x in ORF_list:
    ORF_seq = '' #sequence of full or partial ORF that does not overlap annotated genes, to be added to ORF_list
    overlapping_toggle = False #swtiches to True if the ORF overlaps at all with an annotated gene
    for y in range(min(x[1],x[2])+3, max(x[1],x[2])-3, 3): #exclude start and stop codons
        if sum(genome_coverage[y:y+3]) == 0:
            ORF_seq += gen[y:y+3]
        else:
            overlapping_toggle = True
    if ORF_seq == '':
        x.append('N/A')
        x.append('Overlap')
    else:
        if x[0] == '+':
            x.append(ORF_seq)
        elif x[0] == '-':
            x.append(revcomp(ORF_seq))
        if overlapping_toggle == True:
            x.append('Overlap')
        else:
            x.append('No Overlap')


for x in ORF_list:
    if x[3] != 'N/A':
        pos1 = 0
        pos2 = 0
        pos3 = 0
        for y in range(0, len(x[3]), 3):
            if x[3][y] == 'G' or x[3][y] == 'C':
                pos1 += 1
            if x[3][y + 1] == 'G' or x[3][y + 1] == 'C':
                pos2 += 1
            if x[3][y + 2] == 'G' or x[3][y + 2] == 'C':
                pos3 += 1
        x.append(pos1)
        x.append(pos2)
        x.append(pos3)
        x.append(len(x[3])/3)
        if number_of_tails == 1:
            x.append(scipy.stats.fisher_exact([[pos2, (len(x[3])/3) - pos2],[pos3, (len(x[3])/3) - pos3]],'less')[1])
        elif number_of_tails == 2:
            x.append(scipy.stats.fisher_exact([[pos2, (len(x[3])/3) - pos2],[pos3, (len(x[3])/3) - pos3]],'two-sided')[1])
        


with open(filepath_output_files + output_filename + '.txt', 'w') as output_file:
    output_file.write('Strand\tStart Codon Coordinate\tStop Codon Coordinate\tNon-Overlapping ORF Sequence\tOverlap Status\t')
    output_file.write('G/C Codon Position 1\tG/C Codon Position 2\tG/C Codon Position 3\tTotal Number of Codons\tp-value\n')
    pos1_total = 0
    pos2_total = 0
    pos3_total = 0
    len_total = 0
    for x in ORF_list:
        for y in x:
            output_file.write(str(y) + '\t')
        output_file.write('\n')
        if x[3] != 'N/A':
            pos1_total += int(x[5])
            pos2_total += int(x[6])
            pos3_total += int(x[7])
            len_total += int(x[8])
    output_file.write('\n' + 'Total G/C Codon Position 1' + '\t' + str(pos1_total))
    output_file.write('\n' + 'Total G/C Codon Position 2' + '\t' + str(pos2_total))
    output_file.write('\n' + 'Total G/C Codon Position 3' + '\t' + str(pos3_total))
    output_file.write('\n' + 'Total Codon Count' + '\t' + str(len_total))
    if number_of_tails == 1:
        output_file.write('\n' + str(scipy.stats.fisher_exact([[pos2_total, len_total - pos2_total],[pos3_total, len_total - pos3_total]],'less')[1]))
    elif number_of_tails ==2:
        output_file.write('\n' + str(scipy.stats.fisher_exact([[pos2_total, len_total - pos2_total],[pos3_total, len_total - pos3_total]],'two-sided')[1]))



