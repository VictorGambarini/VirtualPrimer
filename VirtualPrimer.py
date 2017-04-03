#!/usr/bin/python
# VirtualPrimer version 0.1
import os, argparse, math, pandas
import matplotlib.pyplot as plt
from xml.dom import minidom
from collections import OrderedDict
from Bio import SeqIO

################################
#   User Defined Parameters    #
################################
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='[Required]	Input file in fasta formats.', required=True)
parser.add_argument('-p', '--primer', help='[Required]	Primers file in fasta format.', required=True)
parser.add_argument('-o', '--output', help='[Required]	Output folder.', required=True)
parser.add_argument('-l', '--length', help='[Required]	Percentage of length to filter out. (Default = 5)', required=True)
args = parser.parse_args()
length_greater = (int(args.length)+100)/100
length_smaller = (100-(int(args.length)))/100

################################
#    Initial Configuration     #
################################
# Creates the output folders
print("Configuring the program...")
os.system("mkdir -p "+args.output)
os.system("mkdir -p "+args.output+"/tmp/db")
os.system("cp "+args.input+" "+args.output+"/tmp/db/tempdb.fasta")

################################
#       Format Database        #
################################
print("Formating database...")
os.system("makeblastdb -in "+args.output+"/tmp/db/tempdb.fasta -dbtype nucl > /dev/null") 

################################
#       Align Sequences        #
################################
print("Aligning primers...")
os.system("blastn -db "+args.output+"/tmp/db/tempdb.fasta -query "+args.primer+" -out "+args.output+"/tmp/db/out_blast -outfmt 5  -max_hsps 1 -max_target_seqs 1000000000 -word_size 4")

################################
#  Get Data from Blast Output  #
################################
print("Working with blast output...")
out_blast1 = minidom.parse(args.output+"/tmp/db/out_blast")
blastoutput = out_blast1.getElementsByTagName("BlastOutput")[0]
BlastOutput_iterations = blastoutput.getElementsByTagName("BlastOutput_iterations")[0]
Iteration1  = BlastOutput_iterations.getElementsByTagName("Iteration")[0]
Iteration2  = BlastOutput_iterations.getElementsByTagName("Iteration")[1]
Iteration1_hits = Iteration1.getElementsByTagName("Iteration_hits")[0]
Iteration2_hits = Iteration2.getElementsByTagName("Iteration_hits")[0]
Hits1 = Iteration1_hits.getElementsByTagName("Hit")
Hits2 = Iteration2_hits.getElementsByTagName("Hit")
dict1, dict2, dict3 = OrderedDict(), OrderedDict(), OrderedDict()
for Hit in Hits1:
	dict1[Hit.getElementsByTagName("Hit_def")[0].firstChild.data] = (Hit.getElementsByTagName("Hit_len")[0].firstChild.data), (Hit.getElementsByTagName("Hsp_hit-from")[0].firstChild.data), (Hit.getElementsByTagName("Hsp_hit-to")[0].firstChild.data)
for Hit in Hits2:
	dict2[Hit.getElementsByTagName("Hit_def")[0].firstChild.data] = (Hit.getElementsByTagName("Hsp_hit-from")[0].firstChild.data), (Hit.getElementsByTagName("Hsp_hit-to")[0].firstChild.data)
dict3 = OrderedDict((k, dict1[k] + dict2[k]) for k in dict1 if k in dict2)

################################
#   Generate tsv Output File   #
################################
print("Generating output files...")
with open(args.output+'/VirtualPrimer.out', 'w') as tsv_out:
	dict_headers = {}
	dict_headers2 = {}
	for key, value in dict3.items():		
		string_value = (str(value))
		for c in "',)(":
			string_value = string_value.replace(c, '')
		hit1_from = string_value.split(' ')[1]
		hit1_to = string_value.split(' ')[2]
		hit2_from = string_value.split(' ')[3]
		hit2_to = string_value.split(' ')[4]
		if hit1_from < hit2_from:
			if hit1_from < hit1_to:
				hit1 = hit1_from
			else:
				hit1 = hit1_to
			if hit2_from > hit2_to:
				hit2 = hit2_from
			else:
				hit2 = hit2_to	
			hitlength = math.fabs(int(hit1)-int(hit2))
		else:
			if hit1_from > hit1_to:
				hit1 = hit1_from
			else:
				hit1 = hit1_to 
			if hit2_from < hit2_to:
				hit2 = hit2_from
			else:
				hit2 = hit2_to
			hitlength = math.fabs(int(hit1)-int(hit2))
		tsv_out.write(str(key)+'\t'+string_value.replace(' ', '\t')+'\t'+str(hitlength)+'\n')
		dict_headers[str(key).split(" ")[0]] = hitlength
		dict_headers2[str(key).split(" ")[0]+"hit1"] = hit1
		dict_headers2[str(key).split(" ")[0]+"hit2"] = hit2
	length_average = sum(dict_headers.values())/len(dict_headers)

################################
#  Generate Fasta Output File  #
################################
fasta_out = open(args.output+"/fasta_filtered.fa", 'w')
for seq_record in SeqIO.parse(args.input, "fasta"):
	print(seq_record.id)
	if seq_record.id in dict_headers and dict_headers[seq_record.id] <= length_average*length_greater and dict_headers[seq_record.id] >= length_average*length_smaller:
		print(seq_record.id)
		pos1 = int(dict_headers2[seq_record.id+"hit1"])
		pos2 = int(dict_headers2[seq_record.id+"hit2"])
		fasta_out.write(">"+seq_record.id+"\n")
		fasta_out.write(str(seq_record.seq[pos1:pos2])+"\n")
fasta_out.close()		
'''fasta = open(args.input, "fasta")
fasta_variable = fasta.readlines()
fasta_out = open(args.output+"/fasta_filtered.fa", 'w')
switch = False
for line in fasta_variable:
	if line[0] == '>':
		if str(line)[1:-1] in dict_headers and dict_headers[str(line)[1:-1]] <= length_average*length_greater and dict_headers[str(line)[1:-1]] >= length_average*length_smaller:
			fasta_out.write(line)
			switch = True
		else:
			switch = False
	else:
		if switch == True:
			fasta_out.write(line)
fasta.close()
fasta_out.close()'''

################################
#       Generate Graphs        #
################################
print("Generating Graphs...")
#Generate Figure 1 - Alignment
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
x = [1, 2, 3, 4]
Iteration1_stat = Iteration1.getElementsByTagName("Iteration_stat")[0]
Statistics1 = Iteration1_stat.getElementsByTagName("Statistics")[0]
DB_number_of_sequences = Statistics1.getElementsByTagName("Statistics_db-num")[0].firstChild.data
y = [int(DB_number_of_sequences), len(dict1), len(dict2), len(dict3)]
name = ['Initial #\nof Seq.', 'Seq.\naligned P1', 'Seq.\naligned P2', 'Seq.\naligned both']
rects1 = ax1.bar(x, y, align='center', color='g', tick_label=y)
def autolabel(rects):
    for ii,rect in enumerate(rects):
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, '%s'% (name[ii]),
                ha='center', va='bottom')
autolabel(rects1)
plt.ylabel('Number of Sequences')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
plt.savefig(args.output+'/Alignments.png')

#Generate Figure 2 - Length Distribution
colnames = ['Sequence', 'Seq_length', 'P1_from', 'P1_to', 'P2_from', 'P2_to', 'Amplicon_length']
data = pandas.read_csv(args.output+'/VirtualPrimer.out', names=colnames, sep='\t')
Amplicon_lengths = data.Amplicon_length.tolist()
x=(range(0, len(Amplicon_lengths)))
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.scatter(x, Amplicon_lengths)
plt.xlabel('Sequence Number')
plt.ylabel('Length (bp)')
plt.xlim([0,len(Amplicon_lengths)])
plt.savefig(args.output+'/Length_distribution.png')

#Generate Figure 3 - Length Distribution After Filtering
colnames = ['Sequence', 'Seq_length', 'P1_from', 'P1_to', 'P2_from', 'P2_to', 'Amplicon_length']
data = pandas.read_csv(args.output+'/VirtualPrimer.out', names=colnames, sep='\t')
Amplicon_lengths = data.Amplicon_length.tolist()
Amplicon_lengths_final = []
for item in Amplicon_lengths:
	if int(item) <= length_average*length_greater and int(item) >= length_average*length_smaller:
		Amplicon_lengths_final.append(item)
x=(range(0, len(Amplicon_lengths_final)))
fig3 = plt.figure()
ax2 = fig3.add_subplot(111)
ax2.scatter(x, Amplicon_lengths_final)
plt.xlabel('Sequence Number')
plt.ylabel('Length (bp)')
plt.xlim([0,len(Amplicon_lengths_final)])
plt.savefig(args.output+'/Length_distribution2.png')

#Generate Figure 4 - Alignment
fig4 = plt.figure()
ax1 = fig4.add_subplot(111)
x = [1, 2, 3, 4, 5]
Iteration1_stat = Iteration1.getElementsByTagName("Iteration_stat")[0]
Statistics1 = Iteration1_stat.getElementsByTagName("Statistics")[0]
DB_number_of_sequences = Statistics1.getElementsByTagName("Statistics_db-num")[0].firstChild.data
y = [int(DB_number_of_sequences), len(dict1), len(dict2), len(dict3), len(Amplicon_lengths_final)]
name = ['Initial #\nof Seq.', 'Seq.\naligned P1', 'Seq.\naligned P2', 'Seq.\naligned both', 'Seq.\nfiltered']
rects1 = ax1.bar(x, y, align='center', color='g', tick_label=y)
def autolabel(rects):
    for ii,rect in enumerate(rects):
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, '%s'% (name[ii]),
                ha='center', va='bottom')
autolabel(rects1)
plt.ylabel('Number of Sequences')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
plt.savefig(args.output+'/Alignments_filtered.png')


################################
#     Generate HTML Report     #
################################
print("Generating HTML Report...")
os.system("cp img/* "+args.output+"/")
with open(args.output+'/VirtualPrimer_Report.html', 'w') as report_html:
	html = """
	<html>
	<head>
	<title>VirtualPrimer</title>
	</head>
	<body>
	<img src=ViP_Report.png width=20% align=left>
	<br><br><br><br><br><br><br>
	<center><h1><font color="orange">Aligned Sequences</font></h1>
	<img src=Alignments.png width=55%>
	<br><br>
	<center><h1><font color="orange">Length Distribution of Aligned Sequences</font></h1>
	<img src=Length_distribution.png width=55%>
	<br><br>
	<center><h1><font color="orange">Length Distribution Histogram of Aligned Sequences</font></h1>
	<img src=Length_distribution2.png width=55%>
	</center>
	<br><br>
	<center><h1><font color="orange">Filtered Sequences</font></h1>
	<img src=Alignments_filtered.png width=55%>
	</center>
	</body>
	</html>
	"""
	report_html.write(html)
print("Done!")

