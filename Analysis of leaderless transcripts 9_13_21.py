
import matplotlib.pyplot as plt
import pandas
import seaborn
import datetime
import os



date = str(datetime.datetime.now().month) + "_" + str(datetime.datetime.now().day) + "_" + str(datetime.datetime.now().year)
time = str(datetime.datetime.now().hour) + "h_" + str(datetime.datetime.now().minute) + "m"



len_of_genome = 4411532 #length of the Mycobacterium tuberculosis genome





filepath = input("Enter the filepath where you want the output directory to be written\n")
filepath_input_files = input("Enter the filepath where the input files are located\n")

filepath_output_files = filepath + 'Output_files_leaderless_analysis_' + date + '_' + time + '/'
os.makedirs(filepath_output_files)



#Ribo-seq, mapped 3' ends of reads
Ribo_seq_rep1_read_ends_gff = input("Enter filename for Ribo-seq replicate 1 gff\n")
Ribo_seq_rep2_read_ends_gff = input("Enter filename for Ribo-seq replicate 2 gff\n")

#RNA-seq, mapped 3' ends of reads
RNA_seq_rep1_read_ends_gff = input("Enter filename for RNA-seq replicate 1 gff\n")
RNA_seq_rep2_read_ends_gff = input("Enter filename for RNA-seq replicate 2 gff\n")


gene_gff = input("Enter the name of the gene gff file\n")


TSS_file = input("Enter the name of the TSS file\n")


genome_fasta = input("Enter the name of the genome fasta file\n")
f = open(filepath_input_files + genome_fasta, 'r')
qqq = f.readline() #skip the header line
gen = 'x' + f.readline().split()[0].upper() #add an "x" at the beginning so the +1 genome position is position [1]
f.close()


max_length_UTR = int(input("Enter the maximum UTR length to be considered (up to 5 nt)\n"))



##############################################################################################
#                                          FUNCTIONS                                         #
##############################################################################################



#coverage_gff_list is a list of filenames for all coverage gffs to be used
#typically this would be two replicates each for Ribo-seq and RNA-seq
#coverage_gff_names lists labels for the different coverage gff datasets
def get_coverage_score_for_non_overlapping_regions(list_of_ORFs, annotated_ORF_gff_name, coverage_gff_list, len_of_genome, strand_column, start_column, stop_column):

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
        if ORF[0] == '+':
            Toggle = False
            for genepos in range(ORF[start_column], ORF[stop_column] + 1):
                if genome_coverage_annotated_ORFs_plus[genepos] == 0:
                    Toggle = True
                    break
            if Toggle == True:
                list_of_ORFs_with_new_left.append([ORF[0], genepos, ORF[stop_column]])
            elif Toggle == False:
                list_of_ORFs_with_new_left.append(['x']) #if there is complete overlap, still need an entry in the list
                
        elif ORF[0] == '-':
            Toggle = False
            for genepos in range(ORF[stop_column], ORF[start_column] + 1):
                if genome_coverage_annotated_ORFs_minus[genepos] == 0:
                    Toggle = True
                    break
            if Toggle == True:
                list_of_ORFs_with_new_left.append([ORF[0], ORF[start_column], genepos])
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
def get_coverage_score_for_ORFs_without_overlap_concern(list_of_ORFs, coverage_gff_list, len_of_genome, strand_column, start_column, stop_column):


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
            if abs(ORF[stop_column] - ORF[start_column]) > 49:
                coverage = 0
                if ORF[strand_column] == '+':
                    for genpos in range(ORF[start_column], ORF[stop_column] + 1):
                        coverage += genome_coverage_plus[genpos]
                elif ORF[strand_column] == '-':
                    for genpos in range(ORF[stop_column], ORF[start_column] + 1):
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

    leg = plt.legend(loc = "upper left", borderpad = 1, fontsize = 60) #for LL ORF metagene plots, legend needs to be on left
    
    #legend line width
    for x in leg.get_lines():
        x.set_linewidth(10)

    #legend font size
    for x in leg.get_texts():
        x.set_fontsize(60)
    
    #set the y-axis scale
    plt.ylim(0, scale_max)
    plt.xlabel(x_axis_label, fontsize = 60, labelpad = 40)
    plt.ylabel(y_axis_label, fontsize = 60, labelpad = 40)
    plt.gcf().set_size_inches(36,24)

    #set the axis label numbers font size
    plt.tick_params(axis='both', which='major', labelsize=48)
    
    metagene_plot.figure.savefig(filepath_output_files + output_filename + '_' + date + '.png')

    return None


##############################################################################################
##############################################################################################






#reverse complement genome sequence
#adding an 'x' at the start means that position 1 is revcompgen[1]
revcompgen = 'x' + revcomp(gen)


list_of_TSSs_with_metadata = [] #every entry in the list has [strand, coord, first 8 nt]

f = open(filepath_input_files + TSS_file, 'r')

for x in f: #read in TSSs to a list of lists; include coordinate, strand, first 8 nt of transcript sequence
    if x.split()[0] == '+':
        first_8_nt = gen[int(x.split()[1]):int(x.split()[1]) + 8]
    elif x.split()[0] == '-':
        first_8_nt = revcompgen[len_of_genome + 1 - int(x.split()[1]):len_of_genome + 1 - int(x.split()[1]) + 8]
    list_of_TSSs_with_metadata.append([x.split()[0], int(x.split()[1]), first_8_nt])
f.close()



for x in list_of_TSSs_with_metadata: #if there's an ATG or GTG in the 8 nt sequence, append the position of the first ATG/GTG to the list entry
    if x[2].count('ATG') > 0 and x[2].count('GTG') == 0:
        x.append(x[2].find('ATG'))
    elif x[2].count('ATG') == 0 and x[2].count('GTG') > 0:
        x.append(x[2].find('GTG'))
    elif x[2].count('ATG') > 0 and x[2].count('GTG') > 0:
        x.append(min(x[2].find('ATG'), x[2].find('GTG')))
    else:
        x.append('N/A') #Entries with "N/A" at list position [3] don't have an ATG or GTG in the first 8 nt (up to a 5 nt leader)


gene_list = [] #read in gene strands and coordinates from a gff
f = open(filepath_input_files + gene_gff, 'r')
for x in f:
    #append strand, start coord, stop coord, gene ID
    if x.split()[6] == '+':
        gene_list.append([x.split()[6], int(x.split()[3]), int(x.split()[4]), x.split()[8]])
    elif x.split()[6] == '-':
        gene_list.append([x.split()[6], int(x.split()[4]), int(x.split()[3]), x.split()[8]])
f.close()

#adds ORF sequence for any TSS with potential start codon in first 8 nt
for x in list_of_TSSs_with_metadata:
    if x[3] != 'N/A':
        counter = x[3] #when start codon is not at +1, start at start codon

        if x[0] == '+':
            while True: #translate out the sequence from the genome using the TSS coordinate + distance to first ATG/GTG as the starting point
                current_codon = gen[x[1] + counter : x[1] + counter + 3]
                counter += 3
                if current_codon == 'TGA' or current_codon == 'TAA' or current_codon == 'TAG':
                    break
            nt_sequence = gen[x[1] + x[3] : x[1] + counter]
            start_coord = x[1] + x[3]
            stop_coord = x[1] + counter - 1 #substract 1 because last position of stop codon is n+2 not n+3

        elif x[0] == '-':
            while True:
                current_codon = revcompgen[(len_of_genome + 1 - x[1]) + counter : (len_of_genome + 1 - x[1]) + counter + 3]
                counter += 3
                if current_codon == 'TGA' or current_codon == 'TAA' or current_codon == 'TAG':
                    break
            nt_sequence = revcompgen[(len_of_genome + 1 - x[1]) + x[3] : (len_of_genome + 1 - x[1]) + counter]
            start_coord = x[1] - x[3]
            stop_coord = x[1] - counter + 1
        aa_sequence = ORFtranslator(nt_sequence)[0]
        x.append(start_coord)
        x.append(stop_coord)
        x.append(nt_sequence)
        x.append(aa_sequence)



for x in list_of_TSSs_with_metadata: #compare translated out ORFs to annotated ORFs and label as annotated (start codon match), isoform (stop codon match), or novel (remainder)
    if x[3] != 'N/A':
        stop_codon_match_toggle = False
        start_codon_match_toggle = False
        for y in gene_list:
            if y[2] == x[5] and y[0] == x[0]:
                stop_codon_match_toggle = True
            if y[1] == x[4] and y[0] == x[0]:
                start_codon_match_toggle = True
        if start_codon_match_toggle == True:
            x.append('Annotated')
        elif stop_codon_match_toggle == True:
            x.append('Isoform')
        else:
            x.append('Novel')



f = open(filepath_output_files + 'summary of TSSs_' + date + '_' + time + '.txt', 'w') #write the annotated TSS list to a file
f.write('Strand\tTSS coordinate\tFirst 8 nt of transcript\tUTR length (up to 5 nt)\tPotential Start Codon Coordinate\tPotential Stop Codon Coordinate\tnt sequence\taa sequence\tORF classification\n')
for x in list_of_TSSs_with_metadata:
    for y in x:
        f.write(str(y) + '\t')
    f.write('\n')
f.close()



#make metagene lineplot for start codons of non-leaderless TSSs, i.e. TSSs without an ATG/GTG in the first 8 nt
selected_LD_ORFs_for_heatmap = [] #make a shortlist of TSSs that don't have an ATG/GTG in the first 8 nt
for x in list_of_TSSs_with_metadata:
    if x[3] == 'N/A':
        selected_LD_ORFs_for_heatmap.append(x)
lineplot_for_a_single_input_file(selected_LD_ORFs_for_heatmap, 1, 0, 50, 100, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to TSS (nt)', 'Normalized Read Depth',  'leadered_TSSs_two_reps_metagene_plot', 'Greys_r')



#make plots from the lists of potential novel leaderless ORFs
for x in range(max_length_UTR + 1): #change this range to look at start codons at +2, +3, +4, etc.
    selected_ORFs_for_heatmap = []
    for y in list_of_TSSs_with_metadata:
        if y[3] == x:
            selected_ORFs_for_heatmap.append(y)

    selected_novel_ORFs_for_heatmap = []
    for a in selected_ORFs_for_heatmap: #make a sub-list of just the novel ORFs
        toggle = True
        for b in gene_list: #b[2] is annotated gene stop; b[0] is annotated gene strand
            if b[2] == a[5] and b[0] == a[0]: #b[2] == a[5] checks whether a stop codon is in the gene annotation (i.e. annotated or isoform)
                toggle = False
        if toggle == True:
            selected_novel_ORFs_for_heatmap.append(a)

    selected_isoform_ORFs_for_heatmap = []
    for a in selected_ORFs_for_heatmap: #make a sub-list of just the isoform ORFs
        toggle = False
        for b in gene_list:
            if b[1] != a[4] and b[2] == a[5] and b[0] == a[0]: #b[2] == a[5] checks whether a stop codon is in the gene annotation (i.e. annotated or isoform)
                toggle = True
        if toggle == True:
            selected_isoform_ORFs_for_heatmap.append(a)

    selected_annotated_ORFs_for_heatmap = []
    for a in selected_ORFs_for_heatmap: #make a sub-list of just the isoform ORFs
        toggle = False
        for b in gene_list: #b[0] == a[0] matches the strands
            if b[1] == a[4] and b[0] == a[0]: #b[1] == a[4] checks whether a start codon is in the gene annotation (i.e. annotated)
                toggle = True
        if toggle == True:
            selected_annotated_ORFs_for_heatmap.append(a)
    
    if len(selected_novel_ORFs_for_heatmap) > 0:
        lineplot_for_a_single_input_file(selected_novel_ORFs_for_heatmap, 5, 0, 100, 50, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth',  'novel_stop_codons_' + str(x) + '_nt_leader_two_reps_metagene_plot', 'Oranges_r')
        lineplot_for_a_single_input_file(selected_isoform_ORFs_for_heatmap, 5, 0, 100, 50, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth',  'isoform_stop_codons_' + str(x) + '_nt_leader_two_reps_metagene_plot', 'Greens_r')
        lineplot_for_a_single_input_file(selected_annotated_ORFs_for_heatmap, 5, 0, 100, 50, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to Stop Codon (nt)', 'Normalized Read Depth',  'annotated_stop_codons_' + str(x) + '_nt_leader_two_reps_metagene_plot', 'Blues_r')
        lineplot_for_a_single_input_file(selected_novel_ORFs_for_heatmap, 4, 0, 50, 100, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to Start Codon (nt)', 'Normalized Read Depth',  'novel_start_codons_' + str(x) + '_nt_leader_two_reps_metagene_plot', 'Oranges_r')
        lineplot_for_a_single_input_file(selected_isoform_ORFs_for_heatmap, 4, 0, 50, 100, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to Staart Codon (nt)', 'Normalized Read Depth',  'isoform_start_codons_' + str(x) + '_nt_leader_two_reps_metagene_plot', 'Greens_r')
        lineplot_for_a_single_input_file(selected_annotated_ORFs_for_heatmap, 4, 0, 50, 100, 0.2, [filepath_input_files + Ribo_seq_rep1_read_ends_gff, filepath_input_files + Ribo_seq_rep2_read_ends_gff], ['Ribo-seq, Replicate 1', 'Ribo-seq, Replicate 2'], 'Position Relative to Start Codon (nt)', 'Normalized Read Depth',  'annotated_start_codons_' + str(x) + '_nt_leader_two_reps_metagene_plot', 'Blues_r')


#calculate Ribo-seq and RNA-seq coverage across LL novel ORFs, LL annotated ORFs, control set of TSS-->+50 regions
LL_novel_TSS_list_for_coverage_scoring = []
LL_annotated_TSS_list_for_coverage_scoring = []
all_TSS_list_for_coverage_scoring = []



list_of_gffs = []
list_of_gffs.append(filepath_input_files + Ribo_seq_rep1_read_ends_gff)
list_of_gffs.append(filepath_input_files + Ribo_seq_rep2_read_ends_gff)
list_of_gffs.append(filepath_input_files + RNA_seq_rep1_read_ends_gff)
list_of_gffs.append(filepath_input_files + RNA_seq_rep2_read_ends_gff)




for TSS_entry in list_of_TSSs_with_metadata:
    if TSS_entry[3] == 0:
        if TSS_entry[8] == "Novel":
            LL_novel_TSS_list_for_coverage_scoring.append([TSS_entry[0], TSS_entry[4], TSS_entry[5]])
        elif TSS_entry[8] == "Annotated":
            LL_annotated_TSS_list_for_coverage_scoring.append([TSS_entry[0], TSS_entry[4], TSS_entry[5]])


for TSS_entry in list_of_TSSs_with_metadata:
    if TSS_entry[0] == '+':
        all_TSS_list_for_coverage_scoring.append([TSS_entry[0], TSS_entry[1], TSS_entry[1] + 50])
    if TSS_entry[0] == '-':
        all_TSS_list_for_coverage_scoring.append([TSS_entry[0], TSS_entry[1], TSS_entry[1] - 50])



LL_novel_TSS_list_for_coverage_scoring = get_coverage_score_for_non_overlapping_regions(LL_novel_TSS_list_for_coverage_scoring, filepath_input_files + gene_gff, list_of_gffs, len_of_genome, 0, 1, 2)
LL_annotated_TSS_list_for_coverage_scoring = get_coverage_score_for_ORFs_without_overlap_concern(LL_annotated_TSS_list_for_coverage_scoring, list_of_gffs, len_of_genome, 0, 1, 2)

all_TSS_list_for_coverage_scoring = get_coverage_score_for_non_overlapping_regions(all_TSS_list_for_coverage_scoring, filepath_input_files + gene_gff, list_of_gffs, len_of_genome, 0, 1, 2)


f = open(filepath_output_files + 'leaderless_novel_TSS_list_for_coverage_scoring.txt', 'w')
for line in LL_novel_TSS_list_for_coverage_scoring:
    for item in line:
        f.write(str(item) + '\t')
    f.write('\n')
f.close()


f = open(filepath_output_files + 'leaderless_annotated_TSS_list_for_coverage_scoring.txt', 'w')
for line in LL_annotated_TSS_list_for_coverage_scoring:
    for item in line:
        f.write(str(item) + '\t')
    f.write('\n')
f.close()


f = open(filepath_output_files + 'all_TSS_list_for_coverage_scoring.txt', 'w')
for line in all_TSS_list_for_coverage_scoring:
    for item in line:
        f.write(str(item) + '\t')
    f.write('\n')
f.close()

