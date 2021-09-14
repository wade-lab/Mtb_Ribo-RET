

input_gff_file = input("Enter the path and filename for the .gf file to be modified\n")
output_gff_file = input("Enter the path and filename for the output .gff file\n")
mask_file = input("Enter the path and filename for the .txt file listing regions to mask\n")

with open(mask_file, 'r') as exclude:
    newgffplus = []
    newgffminus = []
    a=exclude.readlines()
    for line in a:
        if line.split()[2] == '+':
            newgffplus.append(line.split()[0] + '\t' + line.split()[1] + '\n')
        if line.split()[2] == '-':
            newgffminus.append(line.split()[0] + '\t' + line.split()[1] + '\n')

with open(input_gff_file, 'r') as gff:
    rpmgff = []

    for a in gff:
        toggle = True
        if int(a.split()[5]) > 0:
            for lineplus in newgffplus:
                if int(lineplus.split()[0]) <= int(a.split()[3]) <= int(lineplus.split()[1]):
                    toggle = False
                    break
        if int(a.split()[5]) < 0:
            for lineminus in newgffminus:
                if int(lineminus.split()[0]) <= int(a.split()[3]) <= int(lineminus.split()[1]):
                    toggle = False
                    break
        if toggle == True:
            rpmgff.append(a)

with open(output_gff_file, 'w') as trimmedgff:
    covsum = 0
    for c in rpmgff:
        covsum = covsum + abs(int(c.split()[5]))
    covsum = covsum/1000000

    for c in rpmgff:
        trimmedgff.write(c.split()[0]+'\t'+c.split()[1]+'\tNA3\t'+str(c.split()[3])+'\t'+str(c.split()[4])+'\t'+str(int(c.split()[5])/covsum)+'\t.\t.\t.\n')
 

