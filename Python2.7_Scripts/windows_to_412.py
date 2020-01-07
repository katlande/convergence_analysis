window_rows = []
windows = open("Convert.txt", 'r')
for string in windows:#different VCFs have different end line characters for some reason? This loop checks which the vcf has and parses the lines accordingly
    if "\n" in string:
        window_rows.append(string)
    elif "\r" in string:
        window_rows = string.split("\r")
    elif "\r\n" in string:
        window_rows = string.split("\r\n")
windows.close()
#print window_rows[0]

window_cells = []
for x in window_rows: #For each line in the windows file,
    temp = x.split("\t")#Separate objects in each line into a list where indices start after each \t
    window_cells.append(temp)
#print window_cells[1]#sanity check

full_data_rows = []
fulldata = open("412_gene_loci.txt", 'r')
for string in fulldata:#different VCFs have different end line characters for some reason? This loop checks which the vcf has and parses the lines accordingly
    if "\n" in string:
        full_data_rows.append(string)
    elif "\r" in string:
        full_data_rows = string.split("\r")
    elif "\r\n" in string:
        full_data_rows = string.split("\r\n")
fulldata.close()
#print full_data_rows[1]

full_data_cells = []
for x in full_data_rows: #For each line in the gene info file,
    temp = x.split("\t")#Separate objects in each line into a list where indices start after each \t
    full_data_cells.append(temp)
#print full_data_cells[1]

def single_window_gene(rownum,full_data_cells,windowcells):#this function finds the average rpkm of each tissue within a window. Input rownum (row number of window file), gene nested list (above), window nested list (above)
    matrix = []#empty list to use later
    matrix.append(window_cells[rownum][0])#the first index in the list is the name of the window
    #empty lists to contain tissue counts later:
    isgene = []
    Ha412_name = []
    gene_start_i = []
    gene_end_i = []

    for x in range (1, len(full_data_cells)):#for all the rows of the input gene file,
        
        chr_window = str(window_cells[rownum][1])
        chr_gene = str(full_data_cells[x][3])
        gs = full_data_cells[x][0]
        gene_start = float(gs)
        ge = full_data_cells[x][1]
        gene_end = float(ge)
        ws = window_cells[rownum][2]
        window_start = float(ws)
        we = window_cells[rownum][3]
        window_end = float(we)
        
        #booleans:
        #chromosomes of window and gene match:
        boolean1 = chr_window == chr_gene
        #gene start > window start and < window end:
        boolean2 = gene_start > window_start and gene_start < window_end
        #gene end < window end and gene end > window start:
        boolean3 = gene_start < window_end and gene_start > window_start
        #gene start > window start and gene end < window end:
        boolean4 = gene_start < window_start and gene_end > window_end
        #are any conditions of booleans2,3, or 4 met?
        boolean5 = boolean2 or boolean3 or boolean4
        if boolean1 and boolean5:#so, if a gene overlaps a window, add the gene names to the temporary list for each tissue. List may contain MULTIPLE genes, multiple rpkm values.
            print window_cells[rownum][0]
            Ha412_name.append(full_data_cells[x][2])
            gene_start_i.append(full_data_cells[x][0])
            gene_end_i.append(full_data_cells[x][1])
    
    Ha412_true = []
    window_start = int(window_cells[rownum][2])
    window_end = int(window_cells[rownum][3])
    if len(Ha412_name) > 2:
        for x in range(0,len(Ha412_name)):
            if (int(gene_end_i[x]) > window_start) and (window_end > int(gene_start_i[x])):
                Ha412_true.append(Ha412_name[x])
        isgene.append(len(Ha412_true))
    if len(Ha412_name) == 2:#for each tissue, if the window overlaps with 2 genes:
        if (int(gene_end_i[0]) > window_start) and (window_end > int(gene_start_i[0])) and (int(gene_end_i[1]) > window_start) and (window_end > int(gene_start_i[1])):
            Ha412_true.append(Ha412_name)
            isgene.append(2)
        else:
            if (int(gene_end_i[0]) - window_start) > (window_end - int(gene_start_i[1])):
                Ha412_true.append(Ha412_name[0])
                isgene.append('YES')
            elif (int(gene_end_i[0]) - window_start) < (window_end - int(gene_start_i[1])):
                Ha412_true.append(Ha412_name[1])
                isgene.append('YES')
            elif (int(gene_end_i[0]) - window_start) == (window_end - int(gene_start_i[1])):
                Ha412_true.append(Ha412_name)
                isgene.append('EQUAL COVERAGE')
            else:
                Ha412_true.append('ERROR')
                isgene.append('ERROR')
    elif len(Ha412_name) == 1:
        isgene.append(1)
        Ha412_true.append(Ha412_name)
    elif len(Ha412_name) == 0:#otherwise,
        isgene.append('NA')
        Ha412_true.append('NA')

    matrix.append(isgene)
    matrix.append(Ha412_true)

    return matrix#return the output list

final_results = []
for row in range(1, len(window_cells)):
    final_results.append(single_window_gene(row,full_data_cells,window_cells))

window_data_headers = ['window',"is_gene?", "412_name"]

output = open("windows_to_genes.txt", "w")#add the headers
row_len = len(window_data_headers)
col_len = len(final_results)
for x in range(0, row_len):
    output.write(window_data_headers[x])
    if x < (row_len-1):
        output.write("\t")#spaced by tab
    else:
        output.write("\n")#unless last header
c = 0
while c < col_len:#while loop for row num
    for x in range(0, row_len):#with a nested loop for col num
        s = str(final_results[c][x])
        output.write(s)
        if x < (row_len-1):
            output.write("\t")#same logic - indicies spaced by tabs unless last index in a row
        else:
            output.write("\n")
    c += 1
output.close()