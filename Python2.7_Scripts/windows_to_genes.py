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
fulldata = open("412_to_XRQ_conversion.txt", 'r')
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

#function that averages a list
def Average(lst):
    newlist = []
    for x in lst:
        q = float(x)#make each list object a float
        newlist.append(q)#append it into a new list
    return (sum(newlist)) / (len(newlist))#sum/len of new list so averages have decimals

def single_window_rpkms(rownum,full_data_cells,windowcells):#this function finds the average rpkm of each tissue within a window. Input rownum (row number of window file), gene nested list (above), window nested list (above)
    matrix = []#empty list to use later
    matrix.append(window_cells[rownum][0])#the first index in the list is the name of the window
    #empty lists to contain tissue counts later:
    bractrpkm = []
    corollarpkm = []
    ligulerpkm = []
    ovaryrpkm = []
    pollenrpkm = []
    seedsrpkm = []
    stamenrpkm = []
    stemrpkm = []
    pistilrpkm = []
    isgene = []
    XRQ_name = []
    Ha412_name = []

    for x in range (1, len(full_data_cells)):#for all the rows of the input gene file,
        
        chr_window = str(window_cells[rownum][1])
        chr_gene = str(full_data_cells[x][15])
        gs = full_data_cells[x][13]
        gene_start = float(gs)
        ge = full_data_cells[x][14]
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
        if boolean1 and boolean5:#so, if a gene overlaps a window:
            bractrpkm.append(full_data_cells[x][3])#add the rpkm value to the temporary rpkm list for each tissue. List may contain MULTIPLE genes, multiple rpkm values.
            print window_cells[rownum][0]
            corollarpkm.append(full_data_cells[x][4])
            ligulerpkm.append(full_data_cells[x][5])
            ovaryrpkm.append(full_data_cells[x][6])
            pollenrpkm.append(full_data_cells[x][7])
            seedsrpkm.append(full_data_cells[x][8])
            stamenrpkm.append(full_data_cells[x][9])
            stemrpkm.append(full_data_cells[x][10])
            pistilrpkm.append(full_data_cells[x][11])
            XRQ_name.append(full_data_cells[x][0])
            Ha412_name.append(full_data_cells[x][12])

    if len(bractrpkm) > 0:#for each tissue, if the window overlaps with a gene:
        avgbract = Average(bractrpkm)#average the rpkm value
        isgene.append('YES')#keep track of whether windows contain genes #and which gene it is
    else:#otherwise,
        avgbract = 0.0#set the average rpkm value to 0 (there is no gene here - no gene expression)
        isgene.append('NO')
        XRQ_name.append('NA')
        Ha412_name.append('NA')
    if len(corollarpkm) > 0:
        avgcor= Average(corollarpkm)
    else:
        avgcor = 0.0
    if len(ligulerpkm) > 0:
        avgling= Average(ligulerpkm)
    else:
        avgling = 0.0
    if len(ovaryrpkm) > 0:
        avgov = Average(ovaryrpkm)
    else:
        avgov = 0.0

    if len(pollenrpkm) > 0:
        avgpoll = Average(pollenrpkm)
    else:
        avgpoll = 0.0
    if len(seedsrpkm) > 0:
        avgseed = Average(seedsrpkm)
    else:
        avgseed = 0.0
    if len(stamenrpkm) > 0:
        avgstam = Average(stamenrpkm)
    else:
        avgstam = 0.0
    if len(stemrpkm) > 0:
        avgstem = Average(stemrpkm)
    else:
        avgstem = 0.0
    if len(pistilrpkm) > 0:
        avgpis = Average(pistilrpkm)
    else:
        avgpis = 0.0

    matrix.append(avgbract)#append the average rpkm for the tissue in the window to the output list
    matrix.append(avgcor)
    matrix.append(avgling)
    matrix.append(avgov)
    matrix.append(avgpoll)
    matrix.append(avgseed)
    matrix.append(avgstam)
    matrix.append(avgstem)
    matrix.append(avgpis)
    matrix.append(isgene)
    matrix.append(Ha412_name)
    matrix.append(XRQ_name)

    return matrix#return the output list

final_results = []
for row in range(1, len(window_cells)):
    final_results.append(single_window_rpkms(row,full_data_cells,window_cells))

window_data_headers = ['window', 'Bract.rpkm', 'Corolla.rpkm', 'Ligule.rpkm', 'Ovary.rpkm', 'Pollen.rpkm', 'Seeds.rpkm', 'Stamen.rpkm', 'Stem.rpkm', 'Pistil.rpkm', "is_gene?", "412_name", "XRQ_name"]

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




