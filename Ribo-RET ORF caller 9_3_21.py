
import matplotlib.pyplot as plt
import numpy as np
import pandas
import seaborn
import datetime
from pathlib import Path, PureWindowsPath
import subprocess
import os
import math
import random
import scipy

date = str(datetime.datetime.now().month) + "_" + str(datetime.datetime.now().day) + "_" + str(datetime.datetime.now().year)
time = str(datetime.datetime.now().hour) + "h_" + str(datetime.datetime.now().minute) + "m"



############################
#key parameters and strings#
############################


min_fold = 10 #RPM at potential IERF position must exceed local background by at least "min_fold" times
span = 50 #local background calculated in a window of "span" nt on either side of potential IERF position
threshold = 5.5 #RPM threshold for calling IERFs
window = 10 #potential IERF cannot have a position with higher RPM coverage within a range of "window" nt on either side

#list of acceptable combinations of trinucleotide sequences and positions relative to an IERF
#only these combinations are considered as potential start codons
list_of_good =[['ATG',-14],['ATG',-15],['ATG',-16],['ATG',-17],['ATG',-18],['GTG',-14],['GTG',-15],['GTG',-16],['GTG',-17],['GTG',-18],['TTG',-15]]

len_of_genome = 4411532 #length of the Mtb H37Rv genome


############################
#//////////////////////////#
############################




###########################
#file paths and file names#
###########################



filepath = input("Enter the filepath where you want the output directory to be written\n")
filepath_gffs = input("Enter the filepath where the input files are located\n")

filepath_output_files = filepath + 'Output_files_' + date + '_' + time + '/'
os.makedirs(filepath_output_files)


genome_fasta = input("Enter the name of the genome fasta file\n")

f = open(filepath_gffs + genome_fasta, 'r')
qqq = f.readline() #skip the header line
gen = 'x' + f.readline().split()[0].upper() #add an "x" at the beginning so the +1 genome position is position [1]
f.close()


gene_gff = input("Enter the name of the gene gff file\n")


TSS_file = input("Enter the name of the TSS file\n")



#Ribo-RET, mapped 3' ends of reads
Ribo_RET_rep1_read_ends_gff = input("Enter filename for Ribo-RET replicate 1 gff\n")
Ribo_RET_rep2_read_ends_gff = input("Enter filename for Ribo-RET replicate 2 gff\n")

#Ribo-seq, mapped 3' ends of reads
Ribo_seq_rep1_read_ends_gff = input("Enter filename for Ribo-seq replicate 1 gff\n")
Ribo_seq_rep2_read_ends_gff = input("Enter filename for Ribo-seq replicate 2 gff\n")

#RNA-seq, mapped 3' ends of reads
RNA_seq_rep1_read_ends_gff = input("Enter filename for RNA-seq replicate 1 gff\n")
RNA_seq_rep2_read_ends_gff = input("Enter filename for RNA-seq replicate 2 gff\n")



###########################
#/////////////////////////#
###########################






##############################################################################################
#                                          FUNCTIONS                                         #
##############################################################################################


#coverage_gff_list is a list of filenames for all coverage gffs to be used
#typically this would be two replicates each for Ribo-seq and RNA-seq
#coverage_gff_names lists labels for the different coverage gff datasets
def get_coverage_score_for_non_overlapping_regions(list_of_ORFs, annotated_ORF_gff_name, coverage_gff_list, len_of_genome, strand_col, start_col, stop_col):

    #example ORF entry
    #['+', 'ATG', -15, 3380, 23.48644539, 37.16713881, 'TAG', 3427, 'ATGGTAAGACGAATCTTATTGAGGCACTGTGGTATTCGACGACGTTAG', 'MVRRILLRHCGIRRR*']

    #create lists for each strand, covering all genome positions
    #default value at each position is 0; set to 1 if there is an overlapping ORF with 30 nt extra on either end
    genome_coverage_annotated_ORFs_plus = [0 for i in range(len_of_genome + 51)] # add 1 to count the last nt, since python ranges stop one short of the second number
    genome_coverage_annotated_ORFs_minus = [0 for i in range(len_of_genome + 51)]

    gff_file = open(annotated_ORF_gff_name, 'r')
    for gff_line_unsplit in gff_file:
        gff_line = gff_line_unsplit.split()
        if gff_line[6] == '+':
            left_side = int(gff_line[3])
            right_side = int(gff_line[4])
            for position in range(max(0, left_side - 30), right_side + 31): 
                genome_coverage_annotated_ORFs_plus[position] = 1
        elif gff_line[6] == '-':
            left_side = int(gff_line[3])
            right_side = int(gff_line[4])
            for position in range(max(0, left_side - 30), right_side + 31): 
                genome_coverage_annotated_ORFs_minus[position] = 1
    f.close()



    #find non-overlapping regions
    list_of_ORFs_with_new_left = [] #left-hand side of first part of each ORF that is non-overlapping with an annotated ORF
    
    for ORF in list_of_ORFs:
        if ORF[strand_col] == '+':
            Toggle = False
            for genepos in range(ORF[start_col], ORF[stop_col] + 1):
                if genome_coverage_annotated_ORFs_plus[genepos] == 0:
                    Toggle = True
                    break
            if Toggle == True:
                list_of_ORFs_with_new_left.append([ORF[strand_col], genepos, ORF[stop_col]])
            elif Toggle == False:
                list_of_ORFs_with_new_left.append(['x']) #if there is complete overlap, still need an entry in the list
                
        elif ORF[strand_col] == '-':
            Toggle = False
            for genepos in range(ORF[stop_col], ORF[start_col] + 1):
                if genome_coverage_annotated_ORFs_minus[genepos] == 0:
                    Toggle = True
                    break
            if Toggle == True:
                list_of_ORFs_with_new_left.append([ORF[strand_col], ORF[start_col], genepos])
            elif Toggle == False:
                list_of_ORFs_with_new_left.append(['x'])



    list_of_ORFs_non_overlapping_regions = [] #fix right0hand side to eliminate overlap
    
    for ORF in list_of_ORFs_with_new_left:
        if ORF[0] == '+':
            for genepos in range(ORF[1], ORF[2] + 1):
                if genome_coverage_annotated_ORFs_plus[genepos] == 1:
                    break
            list_of_ORFs_non_overlapping_regions.append([ORF[0], ORF[1], genepos])

        elif ORF[0] == '-':
            for genepos in range(ORF[2], ORF[1] + 1):
                if genome_coverage_annotated_ORFs_minus[genepos] == 1:
                    break
            list_of_ORFs_non_overlapping_regions.append([ORF[0], genepos, ORF[2]])

        else:
            list_of_ORFs_non_overlapping_regions.append(ORF)







    #CALCULATE COVERAGE FOR EACH GFF
    for position_in_list in range(len(coverage_gff_list)):

        #populate coverage scores for each strand
        genome_coverage_plus = [0 for i in range(len_of_genome + 100)]
        genome_coverage_minus = [0 for i in range(len_of_genome + 100)]

        
        gff_file = open(coverage_gff_list[position_in_list], 'r')
        for gff_line in gff_file:
            if float(gff_line.split()[5]) > 0:
                genome_coverage_plus[int(gff_line.split()[3])] = abs(float(gff_line.split()[5]))
            elif float(gff_line.split()[5]) < 0:
                genome_coverage_minus[int(gff_line.split()[3])] = abs(float(gff_line.split()[5]))
        gff_file.close()


        for ORF in list_of_ORFs_non_overlapping_regions:
            Toggle = False
            if len(ORF) > 1:
                if abs(ORF[2] - ORF[1]) > 49:
                    Toggle = True

            if Toggle == True:
                coverage = 0
                if ORF[0] == '+':
                    for genpos in range(ORF[1], ORF[2] + 1):
                        coverage += genome_coverage_plus[genpos]
                elif ORF[0] == '-':
                    for genpos in range(ORF[2], ORF[1] + 1):
                        coverage += genome_coverage_minus[genpos]
                ORF.append(coverage)
            
    for position in range(len(list_of_ORFs)):
        if len(list_of_ORFs_non_overlapping_regions) > 3:
            list_of_ORFs[position].extend(list_of_ORFs_non_overlapping_regions[position][1:])

    return(list_of_ORFs)





##############################################################################################



#coverage_gff_list is a list of filenames for all coverage gffs to be used
#typically this would be two replicates each for Ribo-seq and RNA-seq
#coverage_gff_names lists labels for the different coverage gff datasets
def get_coverage_score_for_ORFs_without_overlap_concern(list_of_ORFs, coverage_gff_list, len_of_genome, strand_col, start_col, stop_col):

    #example ORF entry
    #['+', 'ATG', -15, 3380, 23.48644539, 37.16713881, 'TAG', 3427, 'ATGGTAAGACGAATCTTATTGAGGCACTGTGGTATTCGACGACGTTAG', 'MVRRILLRHCGIRRR*']

    #CALCULATE COVERAGE FOR EACH GFF
    for position_in_list in range(len(coverage_gff_list)):

        #populate coverage scores for each strand
        genome_coverage_plus = [0 for i in range(len_of_genome + 100)]
        genome_coverage_minus = [0 for i in range(len_of_genome + 100)]
        gff_file = open(coverage_gff_list[position_in_list], 'r')
        for gff_line in gff_file:
            if float(gff_line.split()[5]) > 0:
                genome_coverage_plus[int(gff_line.split()[3])] = abs(float(gff_line.split()[5]))
            elif float(gff_line.split()[5]) < 0:
                genome_coverage_minus[int(gff_line.split()[3])] = abs(float(gff_line.split()[5]))
        gff_file.close()

       
        
        for ORF in list_of_ORFs:
            if abs(ORF[stop_col] - ORF[start_col]) > 49:
                coverage = 0
                if ORF[strand_col] == '+':
                    for genpos in range(ORF[start_col], ORF[stop_col] + 1):
                        coverage += genome_coverage_plus[genpos]
                elif ORF[strand_col] == '-':
                    for genpos in range(ORF[stop_col], ORF[start_col] + 1):
                        coverage += genome_coverage_minus[genpos]
                ORF.append(coverage)


    return(list_of_ORFs)





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



def ORFtranslator(sequence):

    aa_codes =  ['*'  ,'*'  ,'*'  ,'I'  ,'I'  ,'I'  ,'L'  ,'L'  ,'L'  ,'L'  ,'L'  ,'L'  ,'V'  ,'V'  ,'V'  ,'V'  ,'F'  ,'F'  ,'M'  ,'C'  ,'C'  ,'A'  ,'A'  ,'A'  ,'A'  ,'G'  ,'G'  ,'G'  ,'G'  ,'P'  ,'P'  ,'P'  ,'P'  ,'T'  ,'T'  ,'T'  ,'T'  ,'S'  ,'S'  ,'S'  ,'S'  ,'S'  ,'S'  ,'Y'  ,'Y'  ,'W'  ,'Q'  ,'Q'  ,'N'  ,'N'  ,'H'  ,'H'  ,'E'  ,'E'  ,'D'  ,'D'  ,'K'  ,'K'  ,'R'  ,'R'  ,'R'  ,'R'  ,'R'  ,'R'  ]
    aa_codons = ['TAA','TAG','TGA','ATT','ATC','ATA','CTT','CTC','CTA','CTG','TTA','TTG','GTT','GTC','GTA','GTG','TTT','TTC','ATG','TGT','TGC','GCT','GCC','GCA','GCG','GGT','GGC','GGA','GGG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','TCT','TCC','TCA','TCG','AGT','AGC','TAT','TAC','TGG','CAA','CAG','AAT','AAC','CAT','CAC','GAA','GAG','GAT','GAC','AAA','AAG','CGT','CGC','CGA','CGG','AGA','AGG']

    amino_acids_dict = {}

    for i in range(0,len(aa_codes)):  
        amino_acids_dict[aa_codons[i]] = aa_codes[i]

    open_reading_frames = []

    for k in range(0,3):
        m = ""
        for j in range(k,len(sequence),3):
            if len(sequence[j:j+3])<3:
                break
            m = m + amino_acids_dict[sequence[j:j+3]]
        open_reading_frames.append(m)
    return(open_reading_frames)

##############################################################################################




def coord_cov_selected(len_of_genome, threshold, window, span, min_fold, filename):
      

    plus_coord = []
    plus_coverage = []

    minus_coord = []
    minus_coverage = []


    plus_strand_coverage = [0 for i in range(len_of_genome +1)]
    minus_strand_coverage = [0 for i in range(len_of_genome +1)]


    plus_coord_selected = []
    plus_coverage_selected = []

    minus_coord_selected = []
    minus_coverage_selected = []

           
    f=open(filename,'r')

    for x in f:
        if float(x.split()[5])<0:
            minus_coord.append(1-(int(x.split()[4]))+(len_of_genome))
            minus_coverage.append(-float(x.split()[5]))
        elif float(x.split()[5])>0:
            plus_coord.append(int(x.split()[4]))
            plus_coverage.append(float(x.split()[5]))
    f.close()


    for x in range(len(plus_coord)):
        plus_strand_coverage[plus_coord[x]] += plus_coverage[x]
    plus_strand_coverage += plus_strand_coverage[0:1000]

    for x in range(len(minus_coord)):
        minus_strand_coverage[minus_coord[x]] += minus_coverage[x]
    minus_strand_coverage += minus_strand_coverage[0:1000]



    for x in range(len(minus_coord)):
        range_values = []
        for y in range(max(minus_coord[x]-span,0),minus_coord[x]+(span+1)):
            range_values.append(minus_strand_coverage[y])
        range_values.sort()
        if minus_coverage[x] < min_fold * (float(sum(range_values))/len(range_values)):
            minus_coverage[x] = 0

    for x in range(len(plus_coord)):
        range_values = []
        for y in range(max(plus_coord[x]-span,0),plus_coord[x]+(span+1)):
            range_values.append(plus_strand_coverage[y])
        range_values.sort()
        if plus_coverage[x] < min_fold * (float(sum(range_values))/len(range_values)):
            plus_coverage[x] = 0
            

    for x in range(len(minus_coord)):
        selection_toggle = False
        if minus_coverage[x] >= threshold:
            selection_toggle = True
        for y in range(minus_coord[x] - window,minus_coord[x] + window + 1):
            if y != minus_coord[x]:
                if minus_strand_coverage[y] > minus_coverage[x]:
                    selection_toggle = False
        if selection_toggle == True:
            minus_coord_selected.append(minus_coord[x])
            minus_coverage_selected.append(minus_coverage[x])

    for x in range(len(plus_coord)):
        selection_toggle = False
        if plus_coverage[x] >= threshold:
            selection_toggle = True
        for y in range(plus_coord[x] - window,plus_coord[x] + window + 1):
            if y != plus_coord[x]:
                if plus_strand_coverage[y] > plus_coverage[x]:
                    selection_toggle = False
        if selection_toggle == True:
            plus_coord_selected.append(plus_coord[x])
            plus_coverage_selected.append(plus_coverage[x])


    return(plus_strand_coverage, plus_coord_selected, plus_coverage_selected, minus_strand_coverage, minus_coord_selected, minus_coverage_selected, filename)

##############################################################################################


def codon_frequency(gen, revcompgen, len_of_genome, shared_peaks_list_plus, shared_peaks_list_minus, codonfilename):


    list_of_sequences_around_peaks = []

    for a in shared_peaks_list_plus:
        if a < 50:
            continue
        list_of_sequences_around_peaks.append(str(gen[a-50:a+10]))

    for b in shared_peaks_list_minus:
        if b < 50:
            continue
        list_of_sequences_around_peaks.append(str(revcompgen[b-50:b+10]))



    nucleotides = ['A','C','G','T']
    codons = []
    for pos1 in nucleotides:
        for pos2 in nucleotides:
            for pos3 in nucleotides:
                codons.append(pos1+pos2+pos3)


    list_of_codon_counts = []
 
    for position in range(len(list_of_sequences_around_peaks[0])-2):
        codon_count = [0 for i in range(len(codons))]
        for t in list_of_sequences_around_peaks:
            trinucleotide = t[position: position+3]
            codon_count[codons.index(trinucleotide)] += 1
        list_of_codon_counts.append(codon_count)


    
    with open(codonfilename + date + '.txt', 'w') as grid:
        for h in range(len(codons)):
            grid.write(str(codons[h]))
            for l in range(len(list_of_codon_counts)):
                grid.write('\t' + str(list_of_codon_counts[l][h]))
            grid.write('\n')

    return(list_of_sequences_around_peaks, list_of_codon_counts)
 

##############################################################################################

#provide a list of lists with ORFs and metadata; specify columns for coordinates and strands for metagene plot
#makes and saves a metagene coverage plot with defined distance window and two coverage gff datasets
#intented to show coverage around start and stop codons for no drug and ret data

def metagene_coverage_dataframe(input_list_of_lists, coordinates_column, strands_column, upstream_distance, downstream_distance, coverage_gff_filename):
    coordinates_and_strands = [] #typically either start codon or stop codon coordinates
    
    #read in coverage numbers from gff --> 2 x dictionary (1 for each strand)
    f = open(coverage_gff_filename, 'r')
    gff_coverage_plus = {}
    gff_coverage_minus = {}
    for x in f:
        if float(x.split()[5]) > 0:
            gff_coverage_plus[int(x.split()[3])] = float(x.split()[5])
        elif float(x.split()[5]) < 0:
            gff_coverage_minus[int(x.split()[3])] = abs(float(x.split()[5]))
    f.close


    #read in coordinates to make metagene plots relative to, and associated strands
    #these are from a gff or other tab-delimited txt file
    #relevant columns are input parameters for the function
    for line in input_list_of_lists:
        coordinates_and_strands.append([line[coordinates_column], line[strands_column]])
    
    
    #lists of lists where the rows correspond to coordinates in the input_list_of_lists file
    #and the columns correspond to genome coordinates in the defined window upstream/downstream
    gff_normalized_coverage_list_of_lists = []

    #populate the lists of lists defined immediately above
    for x in coordinates_and_strands:
        #these temporary lists will be non-normalized, single rows of the final lists of lists
        raw_values_single_row_gff = []
        
        #for plus strand regions, the window goes from [position - upstream distance] to [position + downstream distance + 1]
        if x[1] == '+':
            for y in range(x[0] - upstream_distance, x[0] + downstream_distance + 1):
                if y in gff_coverage_plus:
                    raw_values_single_row_gff.append(gff_coverage_plus[y])
                else:
                    raw_values_single_row_gff.append(0)

        
        #for minus strand regions, the window goes backwards from [position + upstream distance] to [position - downstream distance - 1]
        elif x[1] == '-':
            for y in range(x[0] + upstream_distance, x[0] - downstream_distance -1, -1):
                if y in gff_coverage_minus:
                    raw_values_single_row_gff.append(gff_coverage_minus[y])
                else:
                    raw_values_single_row_gff.append(0)

        
        #the single row values will be normalized to the maximum value in the row
        normalized_values_single_row_gff = []

        
        #if all the values are 0, the row is ignored
        if max(raw_values_single_row_gff) > 0:
            for x in raw_values_single_row_gff:
                normalized_values_single_row_gff.append(x / max(raw_values_single_row_gff))
            gff_normalized_coverage_list_of_lists.append(normalized_values_single_row_gff)
    
    #these are single lists that average the column values from the lists of lists
    averaged_normalized_values_single_row_gff = []

    for x in range(len(gff_normalized_coverage_list_of_lists[0])):
        sum_of_values_in_one_column = 0
        for y in range(len(gff_normalized_coverage_list_of_lists)):
            sum_of_values_in_one_column += gff_normalized_coverage_list_of_lists[y][x]
        averaged_normalized_values_single_row_gff.append(sum_of_values_in_one_column / len(gff_normalized_coverage_list_of_lists))

    return (averaged_normalized_values_single_row_gff)




##############################################################################################

#uses one input file of coordinates/strands and multiple gffs
#input_list_of_lists is a list of start coordinates, stop coordinates, and strands
#list_of_coverage_gff_filenames is a list of all the gffs to be used to make metagene plots from

def lineplot_for_a_single_input_file(input_list_of_lists, coordinates_column, strands_column, upstream_distance, downstream_distance, scale_max, list_of_coverage_gff_filenames, list_of_line_labels, x_axis_label, y_axis_label, output_filename, color):

    df_list = [] #collection of dataframes to make plot

    for coverage_gff_filename in list_of_coverage_gff_filenames:
        df_list.append(metagene_coverage_dataframe(input_list_of_lists, coordinates_column, strands_column, upstream_distance, downstream_distance, coverage_gff_filename))

    #create the x-axis scale for the line-plot
    #upstream and downstream distances are input parameters for the function
    range_values = [x for x in range(-upstream_distance, downstream_distance + 1)]
    
    #create a DataFrame using the x-axis values defined immediately above
    df = pandas.DataFrame(index = range_values)

    #add in the data for gff; behaves like a dictionary
    for x in range(len(list_of_line_labels)):
        df[list_of_line_labels[x]] = df_list[x]

    #clear the previous plot; avoids putting multiple plots on top of each other
    plt.clf()
        
    #create a new plot from the DataFrame defined above
    #dashes = False makes the lines solid rather than dashed
    metagene_plot = seaborn.lineplot(data = df, dashes = False, palette = color, linewidth = 10)

    
    #positions the legend optimally relative to where the data are plotted
    #not sure this actually works since it seems to be undone by next commands
    if upstream_distance > downstream_distance:
        leg = plt.legend(loc = "upper left", borderpad = 1, fontsize = 60)
    else:
        leg = plt.legend(loc = "upper right", borderpad = 1, fontsize = 60)
    
    #legend line width
    for x in leg.get_lines():
        x.set_linewidth(10)

    #legend font size
    for x in leg.get_texts():
        x.set_fontsize(60)
    

    #set the y-axis scale and font size
    plt.ylim(0, scale_max)
    plt.xlabel(x_axis_label, fontsize = 60, labelpad = 40)
    plt.ylabel(y_axis_label, fontsize = 60, labelpad = 40)
    plt.gcf().set_size_inches(36,24)

    #set the axis label numbers font size
    plt.tick_params(axis='both', which='major', labelsize=48)
    
    metagene_plot.figure.savefig(filepath_output_files + output_filename + '_' + date + '.png')

    return None


##############################################################################################

#takes final lists of annotated/novel/isoform ORFs
#start_column specifies the column with start codon coordinate
#strand_column specifies the column with strand
#sequence is extracted from -40 to +20 relative to the start codon
#MFE numbers are calculated by RNAfold
#these numbers are fed to the 'stripplot' function

def MFE_prediction_for_stripplot_input(master_list_of_datasets, start_column, strand_column, labels, y_label, log_or_linear, output_filename):

    #write a fasta file with IDs and sequences for each list in master_list_of_datasets

    #to keep track of the fasta files
    list_of_filenames = []
    
    for x in range(len(master_list_of_datasets)):
        f = open(filepath_output_files + labels[x] + '_RBS_sequences_' + date + '.txt', 'w')
        list_of_filenames.append('"' + filepath_output_files + labels[x] + '_RBS_sequences_' + date + '.txt' + '"')
        for y in master_list_of_datasets[x]:
            if y[strand_column] == '+':
                f.write('>' + str(y[start_column]) + '+\n') #unique ID with coordinate and strand
                f.write(gen[y[start_column] - 40 :y[start_column] + 20] + '\n')
            elif y[strand_column] == '-':
                f.write('>' + str(y[start_column]) + '-\n') #unique ID with coordinate and strand
                #add 1 to make the range the same as on the plus strand
                f.write(revcomp(gen[y[start_column] - 20 +1:y[start_column] + 40 + 1]) + '\n')
        f.close()


    #add a set of 500 random sequences to the master_list_of_datasets
    f = open(filepath_output_files + 'random' + '_RBS_sequences_' + date + '.txt', 'w')
    for x in range(250):
        random_coordinate = random.randint(50, len(gen) - 50)
        f.write('>' + str(random_coordinate) + '+\n') #unique ID
        f.write(gen[random_coordinate - 40: random_coordinate + 20] + '\n')

        random_coordinate = random.randint(50, len(gen) - 50)
        f.write('>' + str(random_coordinate) + '-\n') #unique ID
        f.write(revcomp(gen[random_coordinate - 20 + 1: random_coordinate + 40 + 1]) + '\n')
    f.close()

    list_of_filenames.append('"' + filepath_output_files + 'random' + '_RBS_sequences_' + date + '.txt' + '"')
    labels.append('Random')
    
    
    #location of the RNAfold program
    RNAfold_path = '"C:\\Users\\jtw03\\AppData\\Local\\ViennaRNA Package\\RNAfold.exe"'


    #a list of lists where each sublist is a list of dG values
    list_of_list_of_dGs = []
    
    for x in range(len(list_of_filenames)):
        seq_file = list_of_filenames[x]
        out_file_cmd = '--outfile='
        out_file_name = 'RNAfold_out.txt' #temporary, will be deleted at the end
        RNAfold_output_file = filepath_output_files + labels[x] + '_MFE_deltaG_scores.txt'
        
        #run cmd line for RNAfold then to get default user directory where RNAfold deposits out file
        cmd_list = RNAfold_path + ' -i ' + seq_file + ' ' + out_file_cmd + out_file_name #uses default parameters
        subprocess.run(cmd_list, stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True).stdout.decode('utf8')
        direct = subprocess.run(['dir'], stdout=subprocess.PIPE, shell=True).stdout.decode('utf8').split('\n')
        path = PureWindowsPath(direct[3].rstrip()[14:])
        
        dG_list = []
        
        #get dG from RNAfold out file
        with open(str(Path(path)) + '\\' + out_file_name, 'r') as r, open(RNAfold_output_file, 'w') as w:
            for line in r:
                ID = line.rstrip()[1:]
                seq = next(r).rstrip()
                fold_line = next(r).rstrip()
                dG = float(fold_line.split(' (')[1][0:-1].lstrip(' '))
                w.write(str(dG) + '\n')
                dG_list.append(dG)
        
        #delete postscript and temporary RNAout files
        user_files = os.listdir(str(Path(path)))
        ps_files = [f for f in user_files if f.endswith("_ss.ps")]
        for f in ps_files:
            path_to_file = os.path.join(str(Path(path)), f)
            os.remove(path_to_file)
        os.remove(str(Path(path)) + '\\' + out_file_name)
        
        list_of_list_of_dGs.append(dG_list)


    #Mann-Whitney U test descriptions and p-values
    statistical_comparisons = []
    for x in range(len(list_of_list_of_dGs)-1):
        for y in range(x + 1, len(list_of_list_of_dGs)):
            statistical_comparisons.append(['MFE ' + labels[x] + ' vs ' + labels[y], scipy.stats.mannwhitneyu(list_of_list_of_dGs[x], list_of_list_of_dGs[y], True, 'greater')])

    
    x = stripplot(list_of_list_of_dGs, labels, y_label, log_or_linear, output_filename)

    return statistical_comparisons



##############################################################################################

#takes final lists of annotated/novel/isoform ORFs
#column specifies the column with ret peak height
#these numbers are extracted and fed to the 'stripplot' function

def parse_ret_peak_height_values_for_stripplot_input(master_list_of_datasets, column, labels, y_label, log_or_linear, output_filename):

    list_of_datasets = []
    for x in master_list_of_datasets:
        individual_dataset = [] #just one column from a single list in master_list_of_datasets
        for y in x:
            individual_dataset.append(y[column])
        list_of_datasets.append(individual_dataset)


    #Mann-Whitney U test descriptions and p-values
    statistical_comparisons = []
    for x in range(len(list_of_datasets)-1):
        for y in range(x + 1, len(list_of_datasets)):
            if column == 4:
                statistical_comparisons.append(['Ret peak height rep1 ' + labels[x] + ' vs ' + labels[y], scipy.stats.mannwhitneyu(list_of_datasets[x], list_of_datasets[y], True, 'greater')])
            elif column == 5:
                statistical_comparisons.append(['Ret peak height rep2 ' + labels[x] + ' vs ' + labels[y], scipy.stats.mannwhitneyu(list_of_datasets[x], list_of_datasets[y], True, 'greater')])

    
    x = stripplot(list_of_datasets, labels, y_label, log_or_linear, output_filename)

    return statistical_comparisons


##############################################################################################

#makes a stripplot
#list_of_datasets is a list of lists where each sublist is a set of numbers for one column of the stripplot
#labels is a list of labels that is the same length as list_of_datasets; each entry is the label for a dataset
#y_label is the y-axis label
#log_or_linear sets the y-axis scale to be log or linear

def stripplot(list_of_datasets, labels, y_label, log_or_linear, output_filename):


    combined_data = []
    combined_labels = []
    for x in range(len(list_of_datasets)):
        combined_data += list_of_datasets[x]
        for y in list_of_datasets[x]:
            combined_labels.append(labels[x])

    df = pandas.DataFrame()
    df[y_label] = combined_data
    df[' '] = combined_labels
    
    plt.clf()
    custom_stripplot = seaborn.stripplot(data = df, x = ' ', y = y_label, jitter = 0.4, size = 10) #"size" changes datapoint size
    plt.yscale(log_or_linear)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    
    plt.xlabel(xlabel = '', fontsize = 75, labelpad = 40)
    plt.ylabel(ylabel = y_label, fontsize = 75, labelpad = 40)
    plt.gcf().set_size_inches(24,36)
    plt.gcf().subplots_adjust(bottom = 0.25)
    plt.gcf().subplots_adjust(left = 0.25)

    #set the axis label numbers font size
    plt.tick_params(axis='both', which='major', labelsize=75)
    plt.xticks(rotation = 45)
    
    custom_stripplot.figure.savefig(filepath_output_files + output_filename + date + '.png')
    
    return None


##############################################################################################
#////////////////////////////////////////////////////////////////////////////////////////////#
##############################################################################################

##################
#END OF FUNCTIONS#
##################



#reverse complement genome sequence
#adding an 'x' at the start means that position 1 is revcompgen[1]
revcompgen = 'x' + revcomp(gen)


#import gene position information (start coordinate, stop coordinate, strand) from a gff
#make plus strand start coordinate, plus strand stop coordinate, minus strand start coordinate, and minus strand stop coordinate lists
with open(filepath_gffs + gene_gff, 'r') as genomegff:
    genome_start_coord_plus = []
    genome_start_coord_minus = []
    genome_stop_coord_plus = []
    genome_stop_coord_minus = []
    imported_gene_gff = []
    for gene in genomegff:
        if gene.split()[6] == '+':
            genome_start_coord_plus.append(int(gene.split()[3]))
            genome_stop_coord_plus.append(int(gene.split()[4]))
            imported_gene_gff.append([int(gene.split()[3]), int(gene.split()[4]), gene.split()[6]])
        if gene.split()[6] == '-':
            genome_start_coord_minus.append(int(gene.split()[4]))
            genome_stop_coord_minus.append(int(gene.split()[3]))
            imported_gene_gff.append([int(gene.split()[4]), int(gene.split()[3]), gene.split()[6]])


#call peaks on two replicates of Ribo-RET data, imported as gffs
plus_strand_coverage_retR1, plus_coord_selected_retR1, plus_coverage_selected_retR1, minus_strand_coverage_retR1, minus_coord_selected_retR1, minus_coverage_selected_retR1, filename_retR1 = coord_cov_selected(len_of_genome, threshold, window, span, min_fold, filepath_gffs + Ribo_RET_rep1_read_ends_gff)
plus_strand_coverage_retR2, plus_coord_selected_retR2, plus_coverage_selected_retR2, minus_strand_coverage_retR2, minus_coord_selected_retR2, minus_coverage_selected_retR2, filename_retR2 = coord_cov_selected(len_of_genome, threshold, window, span, min_fold, filepath_gffs + Ribo_RET_rep2_read_ends_gff)


#make lists of peaks common to both replicate Ribo-RET datasets, with a plus strand peak list and a minus strand peak list
shared_peaks_list_plus = []
plus_coord_selected_with_coverage_both_reps = [] #list of lists; each sublist has coordinate, rep1 score, rep2 score, for a shared plus strand peak
for position1 in range(len(plus_coord_selected_retR1)):
    list_to_add = []
  
    for position2 in range(len(plus_coord_selected_retR2)):
        if plus_coord_selected_retR1[position1] == plus_coord_selected_retR2[position2]:
            list_to_add.append(plus_coord_selected_retR1[position1])
            shared_peaks_list_plus.append(plus_coord_selected_retR1[position1])
            list_to_add.append(plus_coverage_selected_retR1[position1])
            list_to_add.append(plus_coverage_selected_retR2[position2])
            plus_coord_selected_with_coverage_both_reps.append(list_to_add)
              
shared_peaks_list_minus = []
minus_coord_selected_with_coverage_both_reps = [] #list of lists; each sublist has coordinate, rep1 score, rep2 score, for a shared minus strand peak
for position1 in range(len(minus_coord_selected_retR1)):
    list_to_add = []
    
    for position2 in range(len(minus_coord_selected_retR2)):
        if minus_coord_selected_retR1[position1] == minus_coord_selected_retR2[position2]:
            list_to_add.append(minus_coord_selected_retR1[position1])
            shared_peaks_list_minus.append(minus_coord_selected_retR1[position1])
            list_to_add.append(minus_coverage_selected_retR1[position1])
            list_to_add.append(minus_coverage_selected_retR2[position2])
            minus_coord_selected_with_coverage_both_reps.append(list_to_add)


start_codon = ['ATG', 'GTG', 'TTG']
coord_position = [-18, -17, -16, -15, -14]
start_codon_coord_plus_shared = []
start_codon_coord_minus_shared = []

for k in plus_coord_selected_with_coverage_both_reps:
    for j in start_codon:
        for i in coord_position:
            if gen[k[0] + i:k[0] + i + 3] == j:
                new_k = [j] + [i] + [k[0] + i] + k[1:len(k)] #adds start codon sequence and coordinate to the sublist to make "new_k"
                start_codon_coord_plus_shared.append(new_k) #make a new list with just the peaks that match start codons

for k in minus_coord_selected_with_coverage_both_reps:
      for j in start_codon:
        for i in coord_position:
            if revcompgen[k[0] + i:k[0] + i + 3] == j:
                new_k = [j] + [i] + [len_of_genome + 1 - k[0] - i] + k[1:len(k)]
                start_codon_coord_minus_shared.append(new_k) #make a new list with just the peaks that match start codons


stop_codon = ['TAA', 'TAG', 'TGA']
revcomp_stop_codon = ['TTA', 'CTA', 'TCA']
revcomp_start_codon = ['CAT', 'CAC', 'CAA']
stop_codon_coord = []
ORF_sequences = []
ORF_sequencesrevcomp = []

gen_extended = gen + gen[1:1000]
revcompgen_extended = revcompgen + revcompgen[1:1000]

for i in range(0,len(start_codon_coord_plus_shared)):
    start = start_codon_coord_plus_shared[i][2]
    for j in range(start, len_of_genome+1000,3):
        trinucleotide = gen_extended[j:j+3]
        if trinucleotide in stop_codon:
            start_codon_coord_plus_shared[i].append(trinucleotide) #appends stop codon sequence
            if j > len_of_genome:
                start_codon_coord_plus_shared[i].append(j+2-len_of_genome) #appends stop codon coordinate, accounting for genome looping
            else:
                start_codon_coord_plus_shared[i].append(j+2) #appends stop codon coordinate
            start_codon_coord_plus_shared[i].append(gen_extended[start:j+3]) #appends nt sequence of ORF
            start_codon_coord_plus_shared[i].append(ORFtranslator(gen_extended[start:j+3])[0]) #appends aa sequence of protein
            break

for i in range(0,len(start_codon_coord_minus_shared)):
    start = start_codon_coord_minus_shared[i][2]
    revcompstart = (len_of_genome + 1) - start
    for m in range(revcompstart, len_of_genome+1000,3):
        trinucleotide = revcompgen_extended[m:m+3]
        if trinucleotide in stop_codon:
            start_codon_coord_minus_shared[i].append(trinucleotide) #appends stop codon sequence
            if m > len_of_genome:
                start_codon_coord_minus_shared[i].append(len_of_genome + 1 - m - 2 - len_of_genome) #appends stop codon coordinate, accounting for genome looping
            else:
                start_codon_coord_minus_shared[i].append(len_of_genome + 1 - m - 2) #appends stop codon coordinate
            start_codon_coord_minus_shared[i].append(revcompgen_extended[revcompstart:m+3]) #appends nt sequence of ORF
            start_codon_coord_minus_shared[i].append(ORFtranslator(revcompgen_extended[revcompstart:m+3])[0]) #appends aa sequence of protein
            break


#make codon frequency table for IERFs
list_of_sequences_around_peaks_shared, lists_of_lists_shared = codon_frequency(gen, revcompgen, len_of_genome, shared_peaks_list_plus, shared_peaks_list_minus, filepath_output_files + 'codon_frequency_shared_peaks')

what_we_want=[]
what_we_dont_want=[]



for ORF in start_codon_coord_plus_shared:
    if ORF[0:2] in list_of_good:
        what_we_want.append(["+"] + ORF)

for ORF in start_codon_coord_minus_shared:
    if ORF[0:2] in list_of_good:
        what_we_want.append(["-"] + ORF)


annotated_ORFs = []
isoform_ORFs = []
novel_ORFs = []

for ORF in what_we_want:
    if ORF[0] == '+':
        if ORF[3] in genome_start_coord_plus:
            annotated_ORFs.append(ORF)
        elif ORF[7] in genome_stop_coord_plus:
            isoform_ORFs.append(ORF)
        else:
            novel_ORFs.append(ORF)
    elif ORF[0] == '-':
        if ORF[3] in genome_start_coord_minus:
            annotated_ORFs.append(ORF)
        elif ORF[7] in genome_stop_coord_minus:
            isoform_ORFs.append(ORF)
        else:
            novel_ORFs.append(ORF)




list_of_gffs = []
list_of_gffs.append(filepath_gffs + Ribo_seq_rep1_read_ends_gff)
list_of_gffs.append(filepath_gffs + Ribo_seq_rep2_read_ends_gff)
list_of_gffs.append(filepath_gffs + RNA_seq_rep1_read_ends_gff)
list_of_gffs.append(filepath_gffs + RNA_seq_rep2_read_ends_gff)

novel_ORFs = get_coverage_score_for_non_overlapping_regions(novel_ORFs, filepath_gffs + gene_gff, list_of_gffs, len_of_genome, 0, 3, 7)
annotated_ORFs = get_coverage_score_for_ORFs_without_overlap_concern(annotated_ORFs, list_of_gffs, len_of_genome, 0, 3, 7)




with open(filepath_output_files + 'Annotated_ORFs_' + date + '.txt', 'w') as annotated_list:
    annotated_list.write('input filename rep 1 = ' + str(filename_retR1)+'\n')
    annotated_list.write('input filename rep 2 = ' + str(filename_retR2)+'\n')
    annotated_list.write('Threshold = ' + str(threshold)+'\n')
    annotated_list.write('Span (local maxima) = ' + str(span)+'\n')
    annotated_list.write('Window = ' + str(window)+'\n')
    annotated_list.write('Min_fold = ' + str(min_fold)+'\n')
    annotated_list.write('Strand\tStart Codon\tStart Codon to Peak Spacing\tStart Coordinate\tRep1 Coverage\tRep2 Coverage\tStop Codon\tStop Coordinate\tORF Sequence\tProtein Sequence\t')
    annotated_list.write('Ribo-seq rep1 coverage\tRibo-seq rep2 coverage\tRNA-seq rep1 coverage\tRNA-seq rep2 coverage\n')

    for ORF in annotated_ORFs:
        for sublist_entry in ORF:
            annotated_list.write(str(sublist_entry) + '\t')
        annotated_list.write('\n')

with open(filepath_output_files + 'Isoform_ORFs_' + date + '.txt', 'w') as isoform_list:
    isoform_list.write('input filename rep 1 = ' + str(filename_retR1)+'\n')
    isoform_list.write('input filename rep 2 = ' + str(filename_retR2)+'\n')
    isoform_list.write('Threshold = ' + str(threshold)+'\n')
    isoform_list.write('Span (local maxima) = ' + str(span)+'\n')
    isoform_list.write('Window = ' + str(window)+'\n')
    isoform_list.write('Min_fold = ' + str(min_fold)+'\n')
    isoform_list.write('Strand\tStart Codon\tStart Codon to Peak Spacing\tStart Coordinate\tRep1 Coverage\tRep2 Coverage\tStop Codon\tStop Coordinate\tORF Sequence\tProtein Sequence\n')

    for ORF in isoform_ORFs:
        for sublist_entry in ORF:
            isoform_list.write(str(sublist_entry) + '\t')
        isoform_list.write('\n')

with open(filepath_output_files + 'Novel_ORFs_' + date + '.txt', 'w') as novel_list:
    novel_list.write('input filename rep 1 = ' + str(filename_retR1)+'\n')
    novel_list.write('input filename rep 2 = ' + str(filename_retR2)+'\n')
    novel_list.write('Threshold = ' + str(threshold)+'\n')
    novel_list.write('Span (local maxima) = ' + str(span)+'\n')
    novel_list.write('Window = ' + str(window)+'\n')
    novel_list.write('Min_fold = ' + str(min_fold)+'\n')
    novel_list.write('Strand\tStart Codon\tStart Codon to Peak Spacing\tStart Coordinate\tRep1 Coverage\tRep2 Coverage\tStop Codon\tStop Coordinate\tORF Sequence\tProtein Sequence\t')
    novel_list.write('Start (overlap removed)\tStop (overlap removed)\tRibo-seq rep1 coverage\tRibo-seq rep2 coverage\tRNA-seq rep1 coverage\tRNA-seq rep2 coverage\n')

    for ORF in novel_ORFs:
        for sublist_entry in ORF:
            novel_list.write(str(sublist_entry) + '\t')
        novel_list.write('\n')



##################################################################
# Make codon frequency tables for annotated, novel, isoform ORFs #
##################################################################

annotated_plus_list = []
annotated_minus_list = []

for x in annotated_ORFs:
    if x[0] == '+':
        annotated_plus_list.append(x[3])
    elif x[0] == '-':
        annotated_minus_list.append(len_of_genome + 1 - x[3])


isoform_plus_list = []
isoform_minus_list = []

for x in isoform_ORFs:
    if x[0] == '+':
        isoform_plus_list.append(x[3])
    elif x[0] == '-':
        isoform_minus_list.append(len_of_genome + 1 - x[3])


novel_plus_list = []
novel_minus_list = []

for x in novel_ORFs:
    if x[0] == '+':
        novel_plus_list.append(x[3])
    elif x[0] == '-':
        novel_minus_list.append(len_of_genome + 1 - x[3])

x,y = codon_frequency(gen, revcompgen, len_of_genome, annotated_plus_list, annotated_minus_list, filepath_output_files + 'annotated_ORF_starts_codon_frequency')
x,y = codon_frequency(gen, revcompgen, len_of_genome, isoform_plus_list, isoform_minus_list, filepath_output_files + 'isoform_ORF_starts_codon_frequency')
x,y = codon_frequency(gen, revcompgen, len_of_genome, novel_plus_list, novel_minus_list, filepath_output_files + 'novel_ORF_starts_codon_frequency')


##################################################################
#////////////////////////////////////////////////////////////////#
##################################################################




####################
# Make Strip Plots #
####################


#make stripplots of ret peak scores for 2x replicates
#statistical comparison is returned
#rep1
Ret_peak_height_rep1_statistical_comaprison = parse_ret_peak_height_values_for_stripplot_input([annotated_ORFs, novel_ORFs, isoform_ORFs], 4, ['Annotated', 'Novel', 'Isoform'], 'Normalized Read Depth (RPM)', 'log', 'Annotated_Novel_Isoform_ret_peak_height_stripplot_rep1')
#rep2
Ret_peak_height_rep2_statistical_comaprison = parse_ret_peak_height_values_for_stripplot_input([annotated_ORFs, novel_ORFs, isoform_ORFs], 5, ['Annotated', 'Novel', 'Isoform'], 'Normalized Read Depth (RPM)', 'log', 'Annotated_Novel_Isoform_ret_peak_height_stripplot_rep2')


#make stripplots of MFE predictions
#statistical comparison is returned
MFE_statistical_comaprison = MFE_prediction_for_stripplot_input([annotated_ORFs, novel_ORFs, isoform_ORFs], 3, 0, ['Annotated', 'Novel', 'Isoform'], 'Delta G', 'linear', 'Annotated_Novel_Isoform_MFE_scores')


####################
#//////////////////#
####################





##########################
# Make list of mock ORFs #
##########################


#control ORFs from RTGs not in any of the called ORF lists
#plus strand only, for simplicity
#rr is a dictionary for genomic positions with a minimum coverage in RNA-seq data
rr = {}
f = open(filepath_gffs + RNA_seq_rep1_read_ends_gff, 'r')
for x in f:
    if float(x.split()[5]) > 0: #ensures plus strand only
        rr[int(x.split()[3])] = float(x.split()[5]) #if there is RNA coverage, it gets put in the dictionary
f.close()


#list of stop codons that match annotated [complete annotation] or novel ORFs
#the mock ORF list should not match any known/predicted or newly discovered ORF at either start or stop codon
disallowed_stops = []
for stop_codon_position in genome_stop_coord_plus:
    disallowed_stops.append(stop_codon_position)

for ORF in novel_ORFs:
    if ORF[0] == '+':
        disallowed_stops.append(ORF[7])


RTG_list_initial = []
counter = 0
while True:
    if counter == 2000: #2000 is enough to ensure that after removing and entries that are present twice (both copies will be deleted)
        break #there will be enough left over to get 1000
    pos = random.randint(1,len_of_genome-1000)
    if pos in rr:
        if gen[pos:pos+3] == 'ATG' or gen[pos:pos+3] == 'GTG':
            if rr[pos] > 0.1: #minimum coverage requirement from RNA-seq data, replicate 1
                ORF_size_counter = 0
                while True:
                    ORF_size_counter += 3
                    current_codon = gen[pos + ORF_size_counter:pos + ORF_size_counter + 3]
                    if current_codon == 'TGA' or current_codon == 'TAA' or current_codon == 'TAG':
                        break
                if pos + ORF_size_counter + 2 not in disallowed_stops:
                    counter += 1
                    RTG_list_initial.append([pos, pos + ORF_size_counter + 2, '+'])


RTG_list_long = [] #longer than it needs to be
RTG_list = [] #the first 1000 from RTG_list_long
for x in RTG_list_initial:
    if RTG_list_initial.count(x) == 1: #throw out duplicates
        RTG_list_long.append(x)

for x in range(1000):
    RTG_list.append(RTG_list_long[x]) #limit the list to 1000


#Write the list of mock ORFs to a file
f = open(filepath_output_files + 'RTG_list.txt', 'w')
f.write('RTG trinucleotide sequences on the + strand \n')
f.write('Starts: ATG, GTG\n')
f.write('Stops: TAA, TGA, TAG\n')
f.write('Start\tStop\tStrand\n')
for x in RTG_list:
    f.write(str(x[0]) + '\t' + str(x[1]) + '\t+\n')
f.close()


##########################
#////////////////////////#
##########################




########################################################
# Make Lists of Leadered and Leaderless Annotated ORFs #
########################################################


#ORFs from the annotation that do not have a TSS within 5 nt
#and also do not have an isoform ORF with a TSS within 5 nt
high_confidence_leadered_annotated_ORFs = []

novel_leaderless_ORFs_0_nt = [] #novel ORFs from TSS analysis with a TSS matching an RTG
annotated_leaderless_ORFs_0_nt = [] #annotated ORFs from TSS analysis with a TSS matching an RTG
isoform_leaderless_ORFs_0_nt = [] #isoform ORFs from TSS analysis with a TSS matching an RTG



potential_leaderless_ORF_stop_codons_plus = [] #list of stop codon positions for ORFs from TSS analysis with a TSS 0-5 nt upstream of an RTG
potential_leaderless_ORF_stop_codons_minus = [] #list of stop codon positions for ORFs from TSS analysis with a TSS 0-5 nt upstream of an RTG


f = open(filepath_gffs + TSS_file, 'r')
for line in f:
    
    if len(line.split()) == 9: #only consider TSSs with a potential ORF nearby
        
        if int(line.split()[3]) <= 5: #only consider TSSs with a potential ORF within 5 nt
            
            if line.split()[0] == '+':
                potential_leaderless_ORF_stop_codons_plus.append(int(line.split()[5]))
            elif line.split()[0] == '-':
                potential_leaderless_ORF_stop_codons_minus.append(int(line.split()[5]))
                
        if line.split()[8] == 'Novel':
            if int(line.split()[3]) == 0:
                novel_leaderless_ORFs_0_nt.append([line.split()[0], int(line.split()[4]), int(line.split()[5])]) #strand, start, stop
                
        elif line.split()[8] == 'Annotated':
            if int(line.split()[3]) == 0:
                annotated_leaderless_ORFs_0_nt.append([line.split()[0], int(line.split()[4]), int(line.split()[5])]) #strand, start, stop
                
        elif line.split()[8] == 'Isoform':
            if int(line.split()[3]) == 0:
                isoform_leaderless_ORFs_0_nt.append([line.split()[0], int(line.split()[4]), int(line.split()[5])]) #strand, start, stop

f.close()


for ORF in imported_gene_gff: #[start, stop, strand]
    if ORF[2] == '+':
        if ORF[1] not in potential_leaderless_ORF_stop_codons_plus:
            high_confidence_leadered_annotated_ORFs.append(ORF)
    elif ORF[2] == '-':
        if ORF[1] not in potential_leaderless_ORF_stop_codons_minus:
            high_confidence_leadered_annotated_ORFs.append(ORF)


########################################################
#//////////////////////////////////////////////////////#
########################################################




###################
# Make Line Plots #
###################


#lists of gffs to make metagene plots from

#two replicates of Ribo-seq
list_of_coverage_gff_filenames_both_reps = []
list_of_coverage_gff_filenames_both_reps.append(filepath_gffs + Ribo_seq_rep1_read_ends_gff)
list_of_coverage_gff_filenames_both_reps.append(filepath_gffs + Ribo_seq_rep2_read_ends_gff)

#two replicates of RNA-seq
list_of_RNA_coverage_gff_filenames_both_reps = []
list_of_RNA_coverage_gff_filenames_both_reps.append(filepath_gffs + RNA_seq_rep1_read_ends_gff)
list_of_RNA_coverage_gff_filenames_both_reps.append(filepath_gffs + RNA_seq_rep2_read_ends_gff)

#one replicate each of Ribo-seq and Ribo-RET
list_of_coverage_gff_filenames_rep1 = []
list_of_coverage_gff_filenames_rep1.append(filepath_gffs + Ribo_seq_rep1_read_ends_gff)
list_of_coverage_gff_filenames_rep1.append(filepath_gffs + Ribo_RET_rep1_read_ends_gff)

#one replicate each of Ribo-seq and Ribo-RET
list_of_coverage_gff_filenames_rep2 = []
list_of_coverage_gff_filenames_rep2.append(filepath_gffs + Ribo_seq_rep2_read_ends_gff)
list_of_coverage_gff_filenames_rep2.append(filepath_gffs + Ribo_RET_rep2_read_ends_gff)


#lists of labels to match the lines on the metagene plots; must be the same length as lists of filenames

#two replicates of Ribo-seq
list_of_labels_both_reps = ['No Drug, Replicate 1', 'No Drug, Replicate 2']

#two replicates of RNA-seq
list_of_labels_RNA_both_reps = ['RNA-seq, Replicate 1', 'RNA-seq, Replicate 2']

#one replicate each of Ribo-seq and Ribo-RET
list_of_labels_rep1 = ['No Drug, Replicate 1', 'Retapamulin, Replicate 1']

#one replicate each of Ribo-seq and Ribo-RET
list_of_labels_rep2 = ['No Drug, Replicate 2', 'Retapamulin, Replicate 2']



#start codon metagene plot for annotated, leadered genes (or any gene set read in from file)
#Ribo-seq vs Ribo-RET
#one replicate per graph
x = lineplot_for_a_single_input_file(high_confidence_leadered_annotated_ORFs, 0, 2, 50, 100, 0.6, list_of_coverage_gff_filenames_rep1, list_of_labels_rep1, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'leadered_annotated_starts_rep1_no_drug_vs_ret', 'deep')
x = lineplot_for_a_single_input_file(high_confidence_leadered_annotated_ORFs, 0, 2, 50, 100, 0.6, list_of_coverage_gff_filenames_rep2, list_of_labels_rep2, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'leadered_annotated_starts_rep2_no_drug_vs_ret', 'deep')


#leadered annotated ORFs start and stop codon metagene plots from Ribo-seq data (two replicates per plot)
x = lineplot_for_a_single_input_file(high_confidence_leadered_annotated_ORFs, 0, 2, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'leadered_annotated_starts_ret_2reps', 'Blues_r')
x = lineplot_for_a_single_input_file(high_confidence_leadered_annotated_ORFs, 1, 2, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'leadered_annotated_stops_ret_2reps', 'Blues_r')


#leaderless annotated ORFs start and stop codon metagene plots from Ribo-seq data (two replicates per plot)
#for 0, 1, 2 and 3 nt 5' UTRs
x = lineplot_for_a_single_input_file(novel_leaderless_ORFs_0_nt, 1, 0, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'leaderless_0_nt_novel_starts_ret_2reps', 'Oranges_r')
x = lineplot_for_a_single_input_file(novel_leaderless_ORFs_0_nt, 2, 0, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'leaderless_0_nt_novel_stops_ret_2reps', 'Oranges_r')

x = lineplot_for_a_single_input_file(annotated_leaderless_ORFs_0_nt, 1, 0, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'leaderless_0_nt_annotated_starts_ret_2reps', 'Blues_r')
x = lineplot_for_a_single_input_file(annotated_leaderless_ORFs_0_nt, 2, 0, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'leaderless_0_nt_annotated_stops_ret_2reps', 'Blues_r')

x = lineplot_for_a_single_input_file(isoform_leaderless_ORFs_0_nt, 1, 0, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'leaderless_0_nt_isoform_starts_ret_2reps', 'Greens_r')
x = lineplot_for_a_single_input_file(isoform_leaderless_ORFs_0_nt, 2, 0, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'leaderless_0_nt_isoform_stops_ret_2reps', 'Greens_r')


#annotated/isoform/novel discovered ORFs start and stop codon metagene plots from Ribo-seq data (two replicates per plot)
x = lineplot_for_a_single_input_file(novel_ORFs, 3, 0, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'novel_starts', 'Oranges_r')
x = lineplot_for_a_single_input_file(novel_ORFs, 7, 0, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'novel_stops', 'Oranges_r')
x = lineplot_for_a_single_input_file(annotated_ORFs, 3, 0, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'annotated_starts', 'Blues_r')
x = lineplot_for_a_single_input_file(annotated_ORFs, 7, 0, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'annotated_stops', 'Blues_r')
x = lineplot_for_a_single_input_file(isoform_ORFs, 3, 0, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'isoform_starts', 'Greens_r')
x = lineplot_for_a_single_input_file(RTG_list, 0, 2, 50, 100, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'RTG_starts', 'Greys_r')
x = lineplot_for_a_single_input_file(RTG_list, 1, 2, 100, 50, 0.3, list_of_coverage_gff_filenames_both_reps, list_of_labels_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'RTG_stops', 'Greys_r')


#annotated/isoform/novel discovered ORFs start and stop codon metagene plots from RNA-seq data (two replicates per plot)
x = lineplot_for_a_single_input_file(novel_ORFs, 3, 0, 50, 100, 0.3, list_of_RNA_coverage_gff_filenames_both_reps, list_of_labels_RNA_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'RNA_novel_starts', 'Oranges_r')
x = lineplot_for_a_single_input_file(novel_ORFs, 7, 0, 100, 50, 0.3, list_of_RNA_coverage_gff_filenames_both_reps, list_of_labels_RNA_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'RNA_novel_stops', 'Oranges_r')
x = lineplot_for_a_single_input_file(annotated_ORFs, 3, 0, 50, 100, 0.3, list_of_RNA_coverage_gff_filenames_both_reps, list_of_labels_RNA_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'RNA_annotated_starts', 'Blues_r')
x = lineplot_for_a_single_input_file(annotated_ORFs, 7, 0, 100, 50, 0.3, list_of_RNA_coverage_gff_filenames_both_reps, list_of_labels_RNA_both_reps, 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth', 'RNA_annotated_stops', 'Blues_r')
x = lineplot_for_a_single_input_file(isoform_ORFs, 3, 0, 50, 100, 0.3, list_of_RNA_coverage_gff_filenames_both_reps, list_of_labels_RNA_both_reps, 'Position Relative to Start Codon (nt)', 'Normalized Read Depth', 'RNA_isoform_starts', 'Greens_r')


##############
#////////////#
##############



#calculate the total number of unique ORFs in each category
ORF_strand_and_start = [] #unique ID for each ORF: strand and start coordinate
for ORF in annotated_ORFs:
    ORF_strand_and_start.append(ORF[0] + str(ORF[3]))
number_of_annotated_ORFs = len(set(ORF_strand_and_start))

ORF_strand_and_start = [] #unique ID for each ORF: strand and start coordinate
for ORF in isoform_ORFs:
    ORF_strand_and_start.append(ORF[0] + str(ORF[3]))
number_of_isoform_ORFs = len(set(ORF_strand_and_start))

ORF_strand_and_start = [] #unique ID for each ORF: strand and start coordinate
for ORF in novel_ORFs:
    ORF_strand_and_start.append(ORF[0] + str(ORF[3]))
number_of_novel_ORFs = len(set(ORF_strand_and_start))


summary_info = []
summary_info.append(['Total # Peaks', len(plus_coord_selected_with_coverage_both_reps) + len(minus_coord_selected_with_coverage_both_reps)])
summary_info.append(['Total # ORFs', number_of_annotated_ORFs + number_of_isoform_ORFs + number_of_novel_ORFs])
summary_info.append(['# Annotated' , number_of_annotated_ORFs])
summary_info.append(['# Isoform' , number_of_isoform_ORFs])
summary_info.append(['# Novel' , number_of_novel_ORFs])







#Estimate the frequency of a false positive matching a start codon


random_plus_coords = []
random_minus_coords = []

sample_size = 100000

for x in range(int(sample_size/2)):
    random_plus_coords.append(random.randint(20,len_of_genome-20))
    random_minus_coords.append(random.randint(20,len_of_genome-20))


start_codons_found = [] #list of lists with each entry being populated with start codon, gap, e.g. ['ATG', -16]
random_plus_coords_selected = [] #list of lists with each entry corresponding to a plus strand random position matching a start codon
random_minus_coords_selected = [] #list of lists with each entry corresponding to a minus strand random position matching a start codon

for k in random_plus_coords:
    for j in start_codon:
        for i in coord_position:
            if gen[k + i:k + i + 3] == j:
                new_k = [j] + [i]
                if new_k in list_of_good:
                    new_k.append(k + i) #new_k = [start codon sequence, gap, position]
                    random_plus_coords_selected.append(new_k)

for k in random_minus_coords:
      for j in start_codon:
        for i in coord_position:
            if revcompgen[k + i:k + i + 3] == j:
                new_k = [j] + [i]
                if new_k in list_of_good:
                    new_k.append(len_of_genome + 1 - k - i)
                    random_minus_coords_selected.append(new_k)





for coord in random_plus_coords_selected:
    start = coord[2]
    for j in range(start, len_of_genome+1000,3):
        trinucleotide = gen_extended[j:j+3]
        if trinucleotide in stop_codon:
            if j > len_of_genome:
                coord.append(j+2-len_of_genome) #appends stop codon coordinate, accounting for genome looping
            else:
                coord.append(j+2) #appends stop codon coordinate
            break

for coord in random_minus_coords_selected:
    start = coord[2]
    revcompstart = (len_of_genome + 1) - start
    for m in range(revcompstart, len_of_genome+1000,3):
        trinucleotide = revcompgen_extended[m:m+3]
        if trinucleotide in stop_codon:
            if m > len_of_genome:
                coord.append(len_of_genome + 1 - m - 2 - len_of_genome) #appends stop codon coordinate, accounting for genome looping
            else:
                coord.append(len_of_genome + 1 - m - 2) #appends stop codon coordinate
            break


random_annotated_count = 0
random_isoform_count = 0
random_novel_count = 0


for x in random_plus_coords_selected:
    if x[2] in genome_start_coord_plus:
        random_annotated_count += 1
    elif x[3] in genome_stop_coord_plus:
        random_isoform_count += 1
    else:
        random_novel_count += 1

for x in random_minus_coords_selected:
    if x[2] in genome_start_coord_minus:
        random_annotated_count += 1
    elif x[3] in genome_stop_coord_minus:
        random_isoform_count += 1
    else:
        random_novel_count += 1



summary_info.append(['# random coordinates with matching start codon', len(random_plus_coords_selected) + len(random_minus_coords_selected)])
summary_info.append(['Sample Size', sample_size])

false_positive_frequency = (len(random_plus_coords_selected) + len(random_minus_coords_selected))/sample_size
N = len(plus_coord_selected_with_coverage_both_reps) + len(minus_coord_selected_with_coverage_both_reps)
total_true_and_false_positives = number_of_annotated_ORFs + number_of_isoform_ORFs + number_of_novel_ORFs
summary_info.append(['FDR', (false_positive_frequency/(1-false_positive_frequency))*(N-total_true_and_false_positives)/total_true_and_false_positives])


fraction_of_false_positives_that_are_annotated = random_annotated_count/(random_annotated_count + random_isoform_count + random_novel_count)
fraction_of_false_positives_that_are_isoform = random_isoform_count/(random_annotated_count + random_isoform_count + random_novel_count)
fraction_of_false_positives_that_are_novel = random_novel_count/(random_annotated_count + random_isoform_count + random_novel_count)

total_number_of_false_positives = (false_positive_frequency/(1-false_positive_frequency))*(N-total_true_and_false_positives)

FDR_annotated = (total_number_of_false_positives * fraction_of_false_positives_that_are_annotated) / number_of_annotated_ORFs
FDR_isoform = (total_number_of_false_positives * fraction_of_false_positives_that_are_isoform) / number_of_isoform_ORFs
FDR_novel = (total_number_of_false_positives * fraction_of_false_positives_that_are_novel) / number_of_novel_ORFs

summary_info.append(['Annotated FDR', FDR_annotated])
summary_info.append(['Isoform FDR', FDR_isoform])
summary_info.append(['Novel FDR', FDR_novel])


#add statistical comparisons
summary_info += MFE_statistical_comaprison
summary_info += Ret_peak_height_rep1_statistical_comaprison
summary_info += Ret_peak_height_rep2_statistical_comaprison

with open(filepath_output_files + 'Run_summary_' + date + "_" + time + '.txt', 'w') as f:
    for x in summary_info:
        for y in x:
            f.write(str(y) + '\t')
        f.write('\n')

with open(filepath_output_files + 'List_of_ERFs_' + date + "_" + time + '.txt', 'w') as f:
    for x in shared_peaks_list_plus:
        f.write(str(x) + '\t+\n')
    for x in shared_peaks_list_minus:
        f.write(str(len_of_genome + 1 - x) + '\t-\n')


