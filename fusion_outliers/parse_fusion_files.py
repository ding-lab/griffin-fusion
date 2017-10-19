import sys

input_files=open(sys.argv[1],'r')
gene_list_file=open(sys.argv[2],'r')

sample_list = []
tool_list = []
fusion_file_list = []
for line in input_files:
  line = line.strip().split()
  sample_list.append(line[0])
  tool_list.append(line[1])
  fusion_file_list.append(line[2])
input_files.close()

gene_list=[]
for gene in gene_list_file:
  gene_list.append(gene.strip())
gene_list_file.close()

for sample_number in range(len(sample_list)):
  sample = sample_list[sample_number]
  tool = tool_list[sample_number]
  fusion_file = open( fusion_file_list[sample_number], 'r')
  fusion_dict = {}
  if tool == "star-fusion":
    fusion_file.readline()
    for line in fusion_file:
      line = line.strip().split('\t')
      FusionName, JunctionReadCount, SpanningFragCount, SpliceType, LeftGene, LeftBreakpoint, RightGene, RightBreakpoint = line[0:8]
      geneA, geneB = FusionName.split("--")
      





output_file=open(sys.argv[3],'w')

output_file.close()
