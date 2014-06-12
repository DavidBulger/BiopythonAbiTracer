#!/usr/bin/python
"""
Copyright (C) 2014 David A Bulger

Reads *.abi Sanger Sequencing files.
Generates mutation *.csv mutation database.
Searches for mutations in those files.
Evaluates *.abi peak height at that location.
Sequence called heterozygous or homozygous for mutation.
"""

import sys
import os
import csv
from Bio import AbiTracer
from Bio.Seq import Seq

# Set working directory
if len(sys.argv) > 1:
    file_location = sys.argv[1]
else:
    file_location = "C:\\Users\\Bulger\\Documents\\GitHub\\biopython\\Tests\\CSV"

os.chdir(file_location)

# Open CSV Storage
# Generate 2D array to store mutation data (currently stores up to 100)
# Header
mutation = [[0 for x in xrange(4)] for x in xrange(101)] 
mutation[0][0] = "Mutation" # Mutation Name
mutation[0][1] = "Sequence" # 20 bp upstream sequence
mutation[0][2] = "WTbp" # Wild-Type Base at Mutation Site
mutation[0][3] = "MTbp" # Mutant Base at Mutation Site

with open('Mutations.csv', 'rb') as csvfile:
    mreader = csv.reader(csvfile)
    j=-1
    for row in mreader:
        j+=1
        mutation[j] = row

# Display Current List of Mutations
print "Do you want to view currently stored mutations?"
StoreTest = raw_input('Enter Y or N:')
Store = 0
if StoreTest is "Y":
    i = 0
    for n in mutation:
        if mutation[i][0] > 0 and mutation[i][0] is not '0':
            print str(i) + "\t" + str(n)
            i+=1
            Store = i
        else: i+=1
    print "Next Available Storage Site:" + str(Store)

# Input new mutation information
print "Do you want to add a new mutation to the search database?"
NewMutantData = raw_input('Enter Y or N:')
if NewMutantData is "Y":
    q = 1
    while q > 0:
        gktest = raw_input('Name of Mutation:')
        SearchID = None
        for index, item in enumerate(mutation):
            if item[0] == gktest:
                print "Mutation already in database"
                SearchID = 1
                break
        if SearchID is 1:
            continue
        gkseq = raw_input('Enter sequence:')
        if len(gkseq) is not 21:
            print "Sequence not 21 base pairs"
            continue
        gkmtest = raw_input('Enter mutant base:')
        if len(gkmtest) is not 1:
            print "Sequence not 1 base pair"
            continue
        print "Mutation Name:" + gktest
        print "Sequence:" + gkseq
        print "Mutant Base:" + gkmtest
        print "Is this information correct?"
        LoopTest = raw_input('Enter Y or N:')
        if LoopTest is "Y":
            mutation[Store][0] = gktest
            mutation[Store][1] = gkseq[0:19]
            mutation[Store][2] = gkseq[20]
            mutation[Store][3] = gkmtest
            q+=1
            ExitTest = raw_input('Any more mutations to enter - Y or N?')
            if ExitTest is "Y":
                Store+=1
            else:
                q = 0

    # Write CSV Storage File
    print "Do you want to save mutation database?"
    StoreData = raw_input('Enter Y or N:')
    if StoreData is "Y":
        with open('Mutations.csv', 'wb') as csvfile:
            mwriter = csv.writer(csvfile)
            mwriter.writerows(mutation)

# Set working directory
if len(sys.argv) > 1:
    file_location = sys.argv[1]
else:
    file_location = "C:\\Users\\Bulger\\Documents\\GitHub\\biopython\\Tests\\Abi"

os.chdir(file_location)

# File Name
files = os.listdir(file_location)
files.sort()

#Definition for Detailed Outputs
#Forward
def print_details(file_name,data1,data2,data3,data4,phd,WTf,lgkf,gkf,fwo,pos):
    gkbp = seq[lgkf:lgkf+21]           #Sequence 20 bp+mutation
    gkF = gkbp[0:19]                    #20bp before mutation
    gkM = gkbp[20]                #Mutation with marker
    i = lgkf+20
    values = data1[pos[i]],data2[pos[i]],data3[pos[i]],data4[pos[i]]
    valueAv = sum(values)/len(values)
    m=-1
    WT = None
    MT = None
    for j in values:
        m+=1
        if j >= valueAv:
            if fwo[m] is mutation[n][2]:
                WT = '1'
            elif fwo[m] is mutation[n][3]:
                MT = '1'
            else:
                continue
    if WT == '1' and MT == None:
        note = "+/+"
    elif WT == '1' and MT == '1':
        note= "m/+"
    elif WT == None and MT == '1':
        note = "m/m"
    elif WT == None and MT == None:
        note = "Mutation not found"                
    print "\t\t"+mutation[n][0]+"  "+note+"\tG is "+str(data1[pos[i]])+  "   \tN2: "+WTf+"\tQS: "+str(phd[i])+"\t"+gkM+" "+gkF
    print "\t\t\t\tA is "+str(data2[pos[i]])+                   "   \tM:  "+mutation[n][3]
    print "\t\t\t\tT is "+str(data3[pos[i]])
    print "\t\t\t\tC is "+str(data4[pos[i]])
    print "------------------------------------------------------------------------------------------------"
    return
    
#Reverse
def print_rdetails(file_name,data1,data2,data3,data4,phd,WTrF,lgkrF,fwo):
    gkbp = seq[lgkrF-1:lgkrF+20]           #Sequence 20 bp+mutation
    gkF = gkbp[1:21]                    #20bp before mutation
    gkM = gkbp[0]                #Mutation with marker
    i = lgkrF-1
    values = data1[pos[i]],data2[pos[i]],data3[pos[i]],data4[pos[i]]
    valueAv = sum(values)/len(values)
    #print valueAv
    m=-1
    WT = None
    MT = None
    for j in values:
        m+=1
        if j >= valueAv*3/4:
            if fwo[m] is WTrF:
                WT = '1'
            elif fwo[m] is gkmr[n]:
                MT = '1'
            else:
                continue
    if WT == '1' and MT == None:
        note = "+/+"
    elif WT == '1' and MT == '1':
        note= "m/+"
    elif WT == None and MT == '1':
        note = "m/m"
    elif WT == None and MT == None:
        note = "Mutation not found"                
    print "\t\t"+mutation[n][0]+"  "+note+"\tG is "+str(data1[pos[i]])+  "   \tN2: "+WTrF+"\tQS: "+str(phd[i])+"\t"+gkM+" "+gkF
    print "\t\t\t\tA is "+str(data2[pos[i]])+                   "   \tM:  "+mutation[n][3]
    print "\t\t\t\tT is "+str(data3[pos[i]])
    print "\t\t\t\tC is "+str(data4[pos[i]])
    print "------------------------------------------------------------------------------------------------"

# Header
print "File Name\tAllele\t  */*\tEvaluation\tKey\tQS\tMutation Site + 20 bp"
print "________________________________________________________________________________________________"

# for loop in order to batch process files
for shortfn in files:
    # Import file
    fn = os.path.join(file_location, shortfn)
    try:
        for record in AbiTracer.AbiIterator(fn):
            file_name = record.name         #File name
            seq = Seq(record.seq)           #Sequence
            rseq = seq.reverse_complement() #Reverse compliment of seq
            phd = list(record.phd)          #Quality score converted to list
            phdr = phd[::-1]                #Reverse phd list
            fwo = record.FWO                #Filter Wheel Order = base order
            data1 = record.DATA1            # fwo[0] = G
            data2 = record.DATA2            # fwo[1] = A
            data3 = record.DATA3            # fwo[2] = T
            data4 = record.DATA4            # fwo[3] = C
            data4r = data4[::-1]
            pos = record.POS                #Base positions in trace files
            posr = pos[::-1]                #Reverse base positions
    except AbiTracer.ABIVersionError, e:
        print 'error', e, fn
		# continue
    print shortfn

    # gk loop in order to batch process all mutations
    n = -1
    for gklist in mutation:
        n+=1 # counter

        # Query sequence variables
        gkf = Seq(gklist[1].upper()) # million mutation project gk allele forward sequence 20 bp before mutation
        WTf = Seq(gklist[2].upper()) # wild type nucleotide at mutation site 

        # Find sequence
        lgkf = seq.find(gkf)  # lgkf = location of gkf in input file str1
        lgkr = rseq.find(gkf) # lgkr = location of gkr in input file str1

        # Display results
        if lgkf >= 0:
            print_details(file_name,data1,data2,data3,data4,phd,WTf,lgkf,gkf,fwo,pos)
            continue
        elif lgkr >= 0:
            gkr = gkn.reversecomplement()
            gkr_ = gkr[n]
            gkrF = gkr_[1:20]
            WTrF = gkr_[0]
            lgkrF = seq.find(gkrF)
            if lgkrF >= 0:
                print_rdetails(file_name,data1,data2,data3,data4,phd,WTrF,lgkrF,fwo)
                continue
            continue
        else:
            continue
    # Space output between files
    print "________________________________________________________________________________________________\n"
