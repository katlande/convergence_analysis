import numpy

#Run mode inputs:

server_y_n = raw_input("Are you running this script on a server? ")
boolean1 = server_y_n is 'yes' or server_y_n is 'y'
boolean2 = server_y_n is 'no' or server_y_n is 'n'
#print boolean1, boolean2#sanity check
if boolean1:
    server = 1
    print "server mode activated"
elif boolean2:
    server = 0
    print "desktop mode activated"
else:
    print "Fatal Error: Unrecognized Input"
    exit()

#add an option to choose between computer time vs. data accuracy:
mem = raw_input("Do you want to run high memory mode? It will take more computing power, but your results will be more accurate. ")
boolean3 = mem is 'yes' or mem is 'y'
boolean4 = mem is 'no' or mem is 'n'
#print boolean1, boolean2#sanity check
if boolean3:
    mem = 1
    print "high memory mode activated"
elif boolean4:
    mem = 0
    print "low memory mode activated"
else:
    print "Fatal Error: Unrecognized Input"
    exit()

#User inputs:
#This script works for both genes and windows (in this case, from GWAS)
#gene_or_window = raw_input("Are you working with genes or windows? (g/w): ")

if server == 0:#if not on server, open dialog box.
    from Tkinter import Tk
    import sys
    #(1): full list of genes/windows with tissue specific RPKM columns
    print "Please choose a file containing the data for the full gene or window set"
    from tkFileDialog import askopenfilename
    Tk().withdraw() # this keeps an extra window from appearing
    full_library = askopenfilename() # show an open file dialog box

    #(2): list of significant genes discovered in analysis
    print "Please choose a file containing the names of your significant genes or windows"
    from tkFileDialog import askopenfilename
    Tk().withdraw() # this keeps an extra window from appearing
    signif_genes = askopenfilename() # show an open file dialog box
elif server == 1:#if on server, raw inputs.
    full_library = raw_input("Please input the name of your library file ")
    signif_genes = raw_input("Please input the name of your significant hits file ")

#(3): desired name of output file
output_file_name = raw_input("Please type a name for your output file (include a file extension): ")

#(4): desired number of permutations
permutation_num = input("Please type the number of permutations you want to run: ")#

header_list = []#headers for the output file.
header_list.append("Organ")
header_list.append("Enrichment_p")
header_list.append("Depletion_p")
header_list.append("null_mean")
header_list.append("null_std")
header_list.append("sample_mean")
header_list.append("sample_std")
#print header_list#sanity check

#the final output -- will add to this in the final step
output_matrix = []
output_matrix.append(header_list)

#set up the input data:
full_library_obj = open(full_library, 'r')#open the library and read each line into python as a list
linelist = []
for string in full_library_obj:#if/else so all end line characters are accepted
    if "\n" in string:
        linelist.append(string)
    elif "\r" in string:
        linelist = string.split("\r")
    elif "\r\n" in string:
        linelist = string.split("\r\n")
    else:
        print "Fatal Error: Endline Character not recognized."#I'm not sure why this would ever happen but I put it in just in case
        exit()
full_library_obj.close()
#print linelist#sanity check

#make a list of lists for each line in the library, where each list is a row and each sublist an index of the row:
library_list = []
for x in linelist: #For each row in the library file...
        temp = x.split("\t")#Separate objects in each line into a list where indices start after each \t
        library_list.append(temp)
#print library_list[1]#sanity check
#make the significant genes/windows into a list:
signif_list = []
InputDS = open(signif_genes, 'r')#open significant genes file
for row in InputDS:#for each row of the DS file
    temp = row.strip()#remove endline character into a temporary object
    signif_list.append(temp)#and append each row list into a bigger list.
InputDS.close()
#print signif_list[1]#sanity check

#Some small functions to be used in the large function. They calculate:
#Function that pulls significant genes from the total list and (1) sums RPKM values and (2) counts the number of genes that are significant (based on locus name).

def scan_library_once(significant_gene_list,library_list,column,signif_index,input_list):#for one significant gene, scan all genes in the library list and look for a match
#significant_gene_list = list of significant genes,library_list = list of lists of all gene data,column = column of the dataframe we're working on,signif_index = row of significant gene list we're working on,input_list = list to append all rpkms into.
    x = 0
    input_list_rpkm = []
    input_list_name = []
    while x < len(library_list):#while loop that scans all rows in the length of the library
        temp1 = library_list[x][0]
        temp2 = significant_gene_list[signif_index]
#        print temp1, temp2
        if temp1 == temp2:#if the gene in the library matches the significant gene of interest
            input_list_rpkm.append(library_list[x][0])#append rpkm value into ouput
            input_list_rpkm.append(library_list[x][column])#append that gene into the output list
        x += 1
    input_list.append(input_list_rpkm)#append each gene/rpkm pair as a nested list
    return input_list#return the name and rpkm value for tthe significant gene
#qq = []
#print scan_library_once(signif_list,library_list,1,1,qq)
#exit()#sanity check

def signif_rpkm(significant_gene_list,library_list,column):#finds the rpkms of ALL the significant genes for a single tissue and appends them into a list
    signif_rows = (len(significant_gene_list))
    counter = 0
    col_rpkm_list = []
    while counter < signif_rows:#for each significant gene, run the scan_library_once function:
        temp = scan_library_once(significant_gene_list,library_list,column,counter,col_rpkm_list)
        counter += 1
    return col_rpkm_list#return a list with the rpkm for your column of interest for all significant genes
#output is a nested list. list[x] = [windowname, rpkm_col]
#print signif_rpkm(signif_list,library_list,1)
#exit()#sanity check

#find all the significant gene/window rpkms and store them in a nested list. On sublist stores ALL rpkms for each tissue:
list_of_sample_rpkm_lists = []
for x in range(1, len(library_list[1])):#for all the columns in the input (tissues)
    temp_rpkm_list = signif_rpkm(signif_list,library_list,x)#make a temporary list with all the gene names, rpkms
    temp_rpkm = []
    for index in range(0, len(temp_rpkm_list)):
        tt = temp_rpkm_list[index][1]#temporary object that separates out only the rpkms for each index of the list
        temp_rpkm.append(tt)#appends the rpkms into a temporary list
    list_of_sample_rpkm_lists.append(temp_rpkm)#then stores each list of significant rpkms/tissue in an output list
#print list_of_sample_rpkm_lists#sanity check

#Set the number of random samples to measure for each permutation of the null distribution. If high memory mode is activated, more will be selected:
if mem == 1:
    if (len(signif_list)*500) < len(library_list):
        number_of_random_samples = (len(signif_list)*100)#number of random samples to collect for the null dsitribution is the number of samples in the significant gene list. If the library is MUCH bigger than the sample list, sample many more samples in the null distribution to shrink the standard deviation of the null distribution a bit. If it's too small, the distribution will be abnormal.
    elif (len(signif_list)*100) < len(library_list) and (len(signif_list)*500) > len(library_list):#add a few layers for library size here
        number_of_random_samples = (len(signif_list)*20)
    else:
        number_of_random_samples = (len(signif_list))
elif mem == 0:
    if (len(signif_list)*500) < len(library_list):
        number_of_random_samples = (len(signif_list)*50)#same as above, but use half as many random samples for low memory mode.
    elif (len(signif_list)*100) < len(library_list) and (len(signif_list)*500) > len(library_list):
        number_of_random_samples = (len(signif_list)*10)
    else:
        number_of_random_samples = (len(signif_list)/2)

#Function that makes a random gene set for a SINGLE permutation, makes both the sample distribution and the null distribution:
def one_permutation_rpkm_list(signif_list,library_list,column,null_sample_number):#null sample number is the number or random samples to take per permutation, as defined above.
    random_rpkm_list = []#list of randomly sampled RPKMs
    while len(random_rpkm_list) < null_sample_number:#until there is an equal number of random genes as the random sample size:
        list_of_used_indicies = []#empty to prevent repetitions
        r = numpy.random.uniform(1, len(library_list))#change the random number each run through the loop
        r = int(r)#set the random number to integer
        if r in list_of_used_indicies:#if the random number has already been used to pull an index
            empty = 0#do nothing -- we don't want repeats
        else:#otherwise
            c = int(column)
            temp = library_list[r][c]#make the rpkm of that index in the library a temporary object
            random_rpkm_list.append(temp)#and append it into the output list
            list_of_used_indicies.append(r)#and the random number into the list of indices not to use again
    #NOW PERMUTATION OF THE SAMPLE DISTRIBUTION:
    true_rpkm_list = list_of_sample_rpkm_lists[(column-1)]#produce the rpkm list for all your significant genes in the column. This is the list that needs to be permuted.
    if (len(true_rpkm_list)) == len(signif_list):
        print "All significant genes/windows matched to library"
    else:
        print true_rpkm_list, "\n", "len f(x) list =", len(true_rpkm_list), "\n\n", signif_list, "len list = ", len(signif_list)
        print "Fatal Error: Some significant genes do not have a match!"
        exit()#sanity check to make sure everything is in working order!
    random_sample_rpkm_list = []
    number_of_random_samples_sdist = (len(true_rpkm_list)/2)
    while len(random_sample_rpkm_list) < number_of_random_samples_sdist:#until there is an equal number of random genes as the random sample size:
        r = numpy.random.uniform(1, len(signif_list))#change the random number each run through the loop
        r = int(r)#set the random number to integer
        if r in list_of_used_indicies:#if the random number has already been used to pull an index
            empty = 0#do nothing -- we don't want repeats
        else:#otherwise
            temp = true_rpkm_list[r]#make the rpkm of that index in the library a temporary object
            random_sample_rpkm_list.append(temp)#and append it into the output list
            list_of_used_indicies.append(r)#and the random number into the list of indices not to
    #return the random permutation and the sample permutation values as one nested list:
    random_and_sample_rpkm_list = []#make an output list
    random_and_sample_rpkm_list.append(random_rpkm_list)#append the null rpkms
    random_and_sample_rpkm_list.append(random_sample_rpkm_list)#append the sample rpkms
    return random_and_sample_rpkm_list#return the randomly sampled rpkms
#print one_permutation_rpkm_list(signif_list,library_list,1,number_of_random_samples)#sanity check
#returns permutation of null at list[0] and permutation of samples at list[1]

#Function that runs the permutation function the number of times that the user initially requested, finds each run's average, and appends it into a list:
def all_random_rpkm_avgs(signif_list,library_list,column):
    z = 0
    permutation_rpkm_list = []
    permutation_sample_rpkm_list = []
    while z < permutation_num:#while the number of loops is less than the specified permutation number
        temp_permutation_list = one_permutation_rpkm_list(signif_list,library_list,column,number_of_random_samples)#make a temporary list with the rpkms from one random sampling
        temp0 = []#two temporary objects, one for the null and one for the sample
        temp1 = []
        for x in range(0, (len(temp_permutation_list[0])-1)):#for all the null rpkms sampled from one permutation...
            temp_0 = float(temp_permutation_list[0][x])
            temp0.append(temp_0)#append the sampled rpkms in
        for x in range(0, (len(temp_permutation_list[1])-1)):#for all the samples from one permutation...
            temp_1 = float(temp_permutation_list[1][x])#append the sampled rpkms in
            temp1.append(temp_1)
        
        temp_avg_random_rpkm = sum(temp0)/len(temp0)#average the rpkms from the random sampling into one value for the null distribution
        temp_avg_random_sample_rpkm = sum(temp1)/len(temp1)#average the rpkms from the random sampling into one value for the sample distribution
        permutation_rpkm_list.append(temp_avg_random_rpkm)#append that value into a list of average rpkms
        permutation_sample_rpkm_list.append(temp_avg_random_sample_rpkm)#and sample average rpkms
        print "permutation number: ", z
        z += 1#continue loop
    permutation_full_null_and_sample = []
    permutation_full_null_and_sample.append(permutation_rpkm_list)
    permutation_full_null_and_sample.append(permutation_sample_rpkm_list)
    return permutation_full_null_and_sample#return all the averages, list[0] = null distribution averages for one column, list[1] = sample distribution averages for one column
#print all_random_rpkm_avgs(signif_list,library_list,1)#sanity check
#exit()


#This is the main function -- uses all the above functions to run a full permutation test on ONE organ (column):
def permutate_one_organ(signif_list,library_list,column):
#first, find average rpkms and std deviations for null and sample distributions:
    list_of_rpkms_from_column = all_random_rpkm_avgs(signif_list,library_list,column)
    temp_random = list_of_rpkms_from_column[0]
    temp_sample = list_of_rpkms_from_column[1]
    mean_random = sum(temp_random)/len(temp_random)
    mean_sample = sum(temp_sample)/len(temp_sample)
    std_random = numpy.std(temp_random)
    std_sample = numpy.std(temp_sample)
#create a counter -- one for all the random average rpkms > rpkm of significant genes, one for random average rpkms < rpkm of significant genes, one for when they are equal (just in case, though unlikely)
    enriched_p_num = 0#empty counters we will add to later
    delpeted_p_num = 0
    NS_p_num = 0
    for x in range (0, (len(temp_random))):#for each index in the average random rpkm list, add 1 to the >, <, or = significant average rpkm counter:
        temp = temp_random[x]
        #print "our random is", temp, "our average is", avg_signif_rpkm#sanity check
        if temp > mean_sample:#add to enrichment, depletion, or NS counters depending on whether sample mean is bigger, smaller, or equal to the null mean respectively.
            delpeted_p_num += 1.0
        elif temp < mean_sample:
            enriched_p_num += 1.0
        elif temp == mean_sample:
            NS_p_num += 1.0
        else:
            print "Error: Non-number detected"
#Function will then calculate probability of enrichment and probability of depletion by:
#p_enrichment = len(depleted_list)/#permutations
#p_depletion = len(enrichment_list)/#permutations
    enrichment_denominator = (permutation_num-NS_p_num)#remove unused values
    p_enrichment = (delpeted_p_num)/(enrichment_denominator)
    p_depletion = (enriched_p_num)/(enrichment_denominator)
#And function will append these values into the output matrix:
    temp_column_list = []#list of the output values from one organ, in the order of the headers (above):
    organ_name = library_list[0][column]
    temp_column_list.append(organ_name)
    temp_column_list.append(p_enrichment)
    temp_column_list.append(p_depletion)
    temp_column_list.append(mean_random)
    temp_column_list.append(std_random)
    temp_column_list.append(mean_sample)
    temp_column_list.append(std_sample)
    print "Permutation of", organ_name, "finished"
    return temp_column_list

#print permutate_one_organ(signif_list,library_list,column)#sanity check
#Now we run this for all organs:
organ_num = len(library_list[0])#set a paramater equal to the number of columns in the rpkm library (1 + organ number, including a window column)
for x in range (1, organ_num):#run the organ permutation function for all columns, from first organ column to last organ column
    temp_organ_vals = permutate_one_organ(signif_list,library_list,x)
    output_matrix.append(temp_organ_vals)#append the output values into a matrix
print output_matrix

#write an output file:
output = open(output_file_name, "w")#add the headers
row_IND = len(output_matrix[0])
for x in range(0, row_IND):
    output.write(header_list[x])
    if x < (row_IND-1):
        output.write("\t")#spaced by tab
    else:
        output.write("\n")#unless last header
c = 1
while c < len(output_matrix):#while loop for row num
    for x in range(0, row_IND):#with a nested loop for col num
        s = str(output_matrix[c][x])
        output.write(s)
        if x < (row_IND-1):
            output.write("\t")#same logic - indicies spaced by tabs unless last index in a row
        else:
            output.write("\n")
    c += 1
output.close()