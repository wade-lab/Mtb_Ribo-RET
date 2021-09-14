
input_fastq_file = input("Enter the path and filename for the .fastq file to be trimmed\n")
output_fastq_file = input("Enter the path and filename for the output .fastq file\n")

f = open(input_fastq_file,'r')
trimmed_fastq = []
while True:
	a=f.readline()
	if a=='':
		break
	b=f.readline()
	c=f.readline()
	d=f.readline()
	if 45 > b.find('CTGTAGGCACC') > 20:
		trimmed_fastq.append(a)
		trimmed_fastq.append(b[0:b.find('CTGTAGGCACC')]+'\n')
		trimmed_fastq.append(c)
		trimmed_fastq.append(d[0:b.find('CTGTAGGCACC')]+'\n')

f.close()
f=open(output_fastq_file,'w')
f.writelines(trimmed_fastq)
f.close()

