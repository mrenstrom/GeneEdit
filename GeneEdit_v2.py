#
# easyGeneEdit.py
#
# Lauren and Mark   3-19-18
# Dnyanada          11-13-18, 11-16-18, 11-27-18
#
# useage: easyGeneEdit.py [fileNameBase] [primer set name]
#
# must run pear previous to running this script to assemble double-ended reads
# <single-ended reads may be placed in stitched_reads sub-directory
# expects fastq file in sub-directory /stitched_reads/fileBaseName.fastq
# 
# This script will filter gene-edit reads for Quality and exact match of primers
# Good Q + primer reads are trimmed <by primer or by optional length>
# Identical reads are combined and marked with an ID
# Combined reads are written to a fasta file
# Needle is run on the fasta file to compute alignments (ml embos)
# Needle results are loaed and the alignment string is matched to the sequence ID
# Aligned sequences are compared to the reference sequence to identify substitution errors 
# Summary gene_editing results are save to /gene_edit sub-directory 
# 
import sys
import math
import re
import csv
import subprocess
from subprocess import call
from collections import OrderedDict
from collections import defaultdict
from collections import namedtuple
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
#
# Named tuples to store CIGAR string token parse info
#
cigToken = namedtuple('cigToken',['cmd','length','rest'])
CigEntry = namedtuple('CigEntry', ['cmd', 'length'], verbose=False)
#--------------------------------------------------------------------------------
#
# parseInsertCigar - given cigar string call cigFindNextToken to break it
# into commands of form ('M', length)
#
#--------------------------------------------------------------------------------
def cigFindNextToken(cig):
    #
    # look for integers then a single (M,I,D)
    #
    r = re.compile("([0-9]+)([a-zA-Z]+)(.+)")
    m = r.match(cig)
    if m:
        return(cigToken(m.group(2),int(m.group(1)),m.group(3)))
    else:
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(cig)
        if m:
            return(cigToken(m.group(2),int(m.group(1)),""))
        else:
            return(cigToken(False,0,""))
#--------------------------------------------------------------------------------
#
# parseInsertCigar - this routine translates a CIGAR string into a list 
# of cigEntry tuples (cmd,length)  
#
#--------------------------------------------------------------------------------
def parseInsertCigar(cig):
    #
    # 186M8I67M 
    #
    # like to return  (M 186)(I 8) (M 67)
    #
    #(186,M) (8,I) (67M)
    retList = []
    #print("Start ",cig)
    while True:
        ct = cigFindNextToken(cig)
        if (ct.cmd == False):
            return retList
        else:
            cigEntry = CigEntry(ct.cmd,ct.length)
            retList.append(cigEntry)

            cig = ct.rest
            if cig == "":
                return retList

    return False
#--------------------------------------------------------------------------------
#
# checkQuality: how many bases are below Q 
# 
#--------------------------------------------------------------------------------
def checkQuality(qualSeq):
  #
  # count how many bases are below low quality
  #
  fail = 0
  i = 0
  for char in qualSeq:
    q = ord(char) - ord('!')
    if (q < 15):
      fail += 1
  return fail
#--------------------------------------------------------------------------------
#
# change to lower case if base not same
#
#--------------------------------------------------------------------------------
def compMod(goodSeq,testSeq):
    tl = [c for c in testSeq]
    i = 0
    for char in goodSeq:
        if char != tl[i]:
            tl[i] = tl[i].lower()
        i += 1
    return ''.join(tl)
#--------------------------------------------------------------------------------
# check for perfect primers, trim
#
# return trimmed seq or ""
#
#
#--------------------------------------------------------------------------------
def checkPrimers(apt,start_primer,end_primer,refseq):
    startPrimerLoc = apt.find(start_primer)
    endPrimerLoc   = apt.find(end_primer)
    
    if startPrimerLoc != 0:
      return ""

    if endPrimerLoc < (len(start_primer) + len(refseq)/2):
      return ""

    return(apt[len(start_primer) : endPrimerLoc])

#------------------------------------------------------------------------------------------------------------
#
# To print the reference, alignment and indels
# TCCCTTCCTCTTTTCTGCTCACACAGGAAGCCCTGGAAGCTGCTTCCTCAGACATGCCGCTGCTGCTACTGCTGCCCCTGCTTCCCTTCCGTGAGTGGCTGTGG
# ||||||||||||*||||||||||||||||||||*||||||**||*|||||||||||||*|||||||||||||||||||||||        ||||*||||||*||
# TCCCTTCCTCTTcTCTGCTCACACAGGAAGCCCgGGAAGCctCTgCCTCAGACATGCCaCTGCTGCTACTGCTGCCCCTGCT--------GTGAaTGGCTGcGG
#
#------------------------------------------------------------------------------------------------------------

def alignment(csvF,txtfile,group,ct,refseq):
    with open(csvF, 'w') as xl:
        xl = csv.writer(xl, delimiter=",", quotechar = '"', quoting=csv.QUOTE_MINIMAL)
        xl.writerow(["ID","Sequence","Count","Frequency","Match/Deletions","Substitution","Insertion"])
        xl.writerow([" ",refseq])
        with open(txtfile,'w') as fh:
            for (k,seq,count,cigar,cgl,subs,insertions) in group:

                # Filter sequences that do not have insertions and deletions
                if cigar.find("D") > -1 or cigar.find("I") > -1:
                    #
                    # Append subs to cigar
                    #
                    fh.write(">{}\n".format(k))
                    fh.write("{}\n".format(count))
                    fh.write(str(float(count)/float(ct))+"\n")

                    freq = (str(float(count)/float(ct)))
                    sub_list = ",".join(map(str, subs))
                    ins_list = "".join(insertions)

                    #print(k + "\t" + seq)
                    #    
                    # To print the reference according to the alignment.
                    #
                    i = 0
                    j = 0
                    for m in cgl:
                        cmdLength = int(m[1])
                        if m[0] == "M":
                          for t in range(0,cmdLength):
                            if seq[i] != refseq[j]:
                              fh.write(refseq[j])
                            else:
                                fh.write(refseq[j])
                            i += 1
                            j += 1
                        elif m[0] == "D":
                            cmdLength = int(m[1])
                            j = j + cmdLength
                            for u in range(0, cmdLength):
                                fh.write(refseq[u])
                        elif m[0] == 'I':
                            space = int(m[1])
                            for s in range(0,space):
                                fh.write(" ")
                            i = i + cmdLength

                    fh.write("\n")  

                    #
                    # To print alignment: Match = |, D = -, Sub = *
                    #
                    i = 0
                    j = 0
                    for m in cgl:
                        cmdLength = int(m[1])
                        if m[0] == "M":
                          for t in range(0,cmdLength):
                            if seq[i] != refseq[j]:
                              fh.write("*")
                            else:
                                fh.write("|")
                            i += 1
                            j += 1
                        elif m[0] == "D":
                            cmdLength = int(m[1])
                            j = j + cmdLength
                            for u in range(0, cmdLength):
                                fh.write(" ")
                        elif m[0] == 'I':
                            space = int(m[1])
                            for s in range(0,space):
                                fh.write(" ")
                            i = i + cmdLength

                    fh.write("\n")  

                    #
                    # To print alignment: D = -, Substitution = lower case
                    #
                    i = 0
                    j = 0
                    al_seq = ""
                    for m in cgl:
                        cmdLength = int(m[1])
                        if m[0] == "M":
                          for t in range(0,cmdLength):
                            if seq[i] != refseq[j]:
                              fh.write(seq[i].lower())
                              al_seq = al_seq + seq[i].lower()
                            else:
                                fh.write(seq[i])
                                al_seq = al_seq + seq[i]
                            i += 1
                            j += 1
                        elif m[0] == "D":
                            cmdLength = int(m[1])
                            j = j + cmdLength
                            for u in range(0, cmdLength):
                                fh.write('-')
                                al_seq = al_seq + "-"
                        elif m[0] == 'I':
                          insertSeq = seq[i:i+cmdLength]
                          fh.write(insertSeq)
                          al_seq = al_seq + insertSeq
                          i = i + cmdLength

                    fh.write("\n")

                    xl.writerow([k,al_seq,count,freq,cigar,sub_list,ins_list])
                                    
                    #
                    # indel cigar
                    #
                    fh.write("INDEL,")
                    first = True
                    for m in cgl:
                      if first:
                        fh.write("{},{}".format(m[0],m[1]))
                      else:
                        fh.write(",{},{}".format(m[0],m[1]))
                      first = False
                    fh.write("\n")
                    #
                    # substitutions
                    #
                    fh.write("SUB,")
                    first = True
                    for s in subs:
                      if first:
                        fh.write(str(s[0]) + "," + s[1] + "," + s[3])
                      else:
                        fh.write("," + str(s[0]) + "," + s[1] + "," + s[3])
                      first = False
                    fh.write("\n")

        return print("Results written in the files!")

#--------------------------------------------------------------------------------
#
# Make fasta files for alignment; Run needle
# 
#--------------------------------------------------------------------------------
def needle(sd,outFile,refFile,samFile):
    #
    # sort in order of sequence count
    #
    ssd = sorted(sd.items(),key=lambda al: al[1],reverse=True)  # smallest first
    #
    # create an ID for each combined seq group, store in preSeqDict based on ID
    #
    preSeqDict = {}
    #
    # create fasta file for alignment
    #
    with open(outFile,'w') as fh:
        id = 0
        for k,v in ssd:
            if v > 1:
              seqID = "ID_{}".format(id)
              preSeqDict[seqID] = (k,v)
              fh.write(">" + seqID + "\n")
              fh.write("{}\n".format(k))
              id += 1
              
    # Call out to needle to generate alignment string - CIGAR        
    print("Step 3: Run Needle")
    # if there is a problem, remember to "ml emboss"
    #
    print(refFile,outFile,samFile)
    p = subprocess.Popen(["needle", 
        "-asequence",refFile,
        "-bsequence",outFile,
        "-gapopen","10.0",
        "-gapextend","0.5",
        "-aformat3","sam",
        "-outfile",samFile], stdout=subprocess.PIPE)
    print(p.communicate())

    return preSeqDict
#-----------------------------------------------------------------------------
#
# Identify substitutions based on alignment
# Create sequence groups by combining information for each seqeunce.
#
#------------------------------------------------------------------------------
def seq_group(samFile,preSeqDict,refseq):
    #
    # Read in sam alignment, based in ID add to seqDict
    #
    seqDict = {}
    with open(samFile,'r') as fh:
      while True:
        l = fh.readline().strip()
        if l == "": break
        l = l.split('\t')
        if l[0][0] == "@": continue
        k = l[0]
        cigar = l[5]
        cgList = parseInsertCigar(cigar)
        #
        # get original sequence and count from seqDict
        # -- can't find key is fatal error --
        s,count = preSeqDict[k]
        #
        # add on cigar string and cgList (tuple list of alignment)
        #   -- sequences with really long cigar probably mean primer bound
        #      to some remote DNA and alignment is bonkers
        #      12 is arbitrary...some seq ave lots of SNPs but usually 3 is normal
        #
        if len(cgList) <= 12:
          seqDict[k] = (s,count,cigar,cgList)
        else:
          print("Bad alignment for sequence")
          #print(s)
          print(cigar)

    #
    # don't accidentally use
    #
    preSeqDict = None
    #
    # Use alignment + sequence + reference sequence to identify substitutions and insertion sequences 
    #
    # sortD = list if dictionary k,v pairs sorted by count
    #
    sortD = sorted(seqDict.items(),key=lambda al: al[1][1],reverse=True)  # smallest first
    #
    # finalGroups will be the final list of sequence groups
    #
    finalGroups =[]
    #
    # for each sequence group do a base-by-base comparison to refseq 
    #
    for k,(seq,count,cigar,cgl) in sortD:
      #print("012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789")
      #for m in cgl:
      #  print(m) 
      #
      # substitution mutations and insertion sequences
      #
      subs = []
      insertions = []
      align_seq = []
      #
      # base by base comparison
      #    i = test sequence index
      #    j = reference sequence index
      #    cgl = list of (command,length) from CIGAR  M = match, I = insert, D = deletion
      #
      i = 0
      j = 0

      for m in cgl:
        cmdLength = int(m[1])
        if m[0] == "M":
          for t in range(0,cmdLength):
            #print("comp {} {}".format(i,j))
            if seq[i] != refseq[j]:
              #print("substitue at ",i,seq[i],j,refseq[j])
              mut = (j+1,refseq[j],i+1,seq[i])  # Increase the position no. by 1 (i and j) because the count in python begins at 0.
              subs.append(mut)
            i += 1
            j += 1
        elif m[0] == "D":
          j = j + cmdLength
        elif m[0] == 'I':
          insertSeq = seq[i:i+cmdLength]
          insertions.append(insertSeq)
          i = i + cmdLength

     #
     # add entry for sequence group with all data with the condition to limit substitutions to 12.
     # 
      if len(subs) <= 10:
        finalGroups.append((k,seq,count,cigar,cgl,subs,insertions))

    return finalGroups
       
#---------------------------------------------------------------------------
#   Main Entry
#
#
# @M03100:356:000000000-C4NGM:1:1101:17289:2061 1:N:0:2
# CCCACACTATCTCAATGCAAATATCTGTCTGAAACGGTCCCTGGCTAAACTCCACCCAT...
# +
# >>>1>11C1D@DD33B3B111B1FGHBGGHDBHHHCAFGGHHHHHHHHHHHHHHHGGGH...
#
#
#
#---------------------------------------------------------------------------
#
# useage: python filterFasta.py file.fastq out.fastq
#
#
if len(sys.argv) != 3:
  print("python filterFastq.py inputFileiBase primerCode")
  exit(0)
script,inputFileBase,primerCode = sys.argv
#
# should be params
#
if primerCode == "BE3":
  start_primer = "GGATCCAAATTTCTGGCTGC"
  end_primer   = "GATCCCAGTAGGAACAACTGC"
  refseq       = "AAGTGCAGGAGTCAGTGACGGTACAGGAGGGTTTGTGCGTCCTCGTGCCCTGCACTTTCTTCCATCCCATACCCTACTACGACAAGAACTCCCCAGTTCATGGTTACTGGTTCCGGGAAGGAGCCATTATATCCAGGGACTCTCCAGTGGCCACAAACAAGCTAGATCAAGAAGTACAGGAGGAGACTCAGGGCAGATTCCGCCTCCTTGGG"
  refFile      = "BE.fa"
elif primerCode == "CD33E1":
  start_primer = "CTGTAGTCCTTCCCCTCCAC"
  end_primer   = "CTGCAAGTGCAGGAGTCAG"
  refseq       = "TCCCTTCCTCTTTTCTGCTCACACAGGAAGCCCTGGAAGCTGCTTCCTCAGACATGCCGCTGCTGCTACTGCTGCCCCTGCTGTGGGCAGGTGAGTGGCTGTGGGGAGAGGGGTTGTCGGGCTGGGCCGAGCTGACCCTCGTTTCCCCACAGGGGCCCTGGCTATGGATCCAAATTTCTGG"
  refFile      = "CD33E1.fa"
  refFile2     = "Pseudogene.fa"
  #pseudo       = "TGCCACTGCTGCTACTGCTGCCCCTGCTGTGGGCAGGTGAATGGCTGCGGGGAGAGGGGTTGTCGGGCTGGGCCGAGCTGACCCTCGTTTCCCCACAGGGGCCCTGGCTATGGATCCAAAAATCCGGCTGCAAGTGCAGGAGTCAGTGACGGTACAGGAGGGTTTGTGCGTCCTCGTGCCCTGCACTTTCTCCCATCCCATACCCTACTACAACAGGAATTTCTCAGTTCATGGTTACTGGTTCCGGGAAGGAGCCATTGTATCCAGGGACTCTCCAGTGGCCACAAACAAGCTAGATCAAGAAGTACAGGAGGAGACTCAGGGCTGATTCCGCCTCCTTGGGGATCCCAGTAAGAACAACTGCTCCCTGAGCATCGTAGACGCCAGGAGGAGGGATAATCGTTCATACTTCTTTCGGATGGAGAGAGGAAGTACCAAACACAGTTACAAATCTCCCCAGCTCTCTGTGCATGTGACAGGTGAGGCACAGGCTTCAGAAGCAGCCACAAGGGAAGGTCAAAGGGACCTCAGGACAGGGCTTGGGATGGGACCCATGTCCTGGAAGGGGGTTGGGAATGAAGCCTGTC"
  pseudo       = "TGCCACTGCTGCTACTGCTGCCCCTGCTGTGGGCAGGTGAATGGCTGCGGGGAGAGGGGTTGTCGGGCTGGGCCGAGCTGACCCTCGTTTCCCCACAGGGGCCCTGGCTATGGATCCAAAAATCCGG"
else:
  print("Unknown primer code")
  exit(0)
#
# open fastq
#
print("Step 1: Creating output files.")
inFile      = "/home/dpande/stitched_reads/" + inputFileBase + ".fastq"
outFile     = "/home/dpande/try/" + inputFileBase + "_ge.fasta"
pseudoFile  = "/home/dpande/try/" + inputFileBase + "_pseudo_ge.fasta"
samFile     = "/home/dpande/try/" + inputFileBase + "_ge.sam"
pSamFile    = "/home/dpande/try/" + inputFileBase + "_pseudo_ge.sam"
mutFile     = "/home/dpande/try/" + inputFileBase + "_ge.txt"
pmutFile    = "/home/dpande/try/" + inputFileBase + "_pseudo_ge.txt"
csvFile     = "/home/dpande/try/" + inputFileBase + "_ge.csv"
pcsvFile    = "/home/dpande/try/" + inputFileBase + "_pseudogene_ge.csv"

ref = open("/home/dpande/" + refFile, 'r')
ref_hd = ref.readline().strip()
ref_seq = ref.readline().strip()

allCount    = 0
failQual    = 0
failPrimers = 0
goodCount   = 0
pseudoCount = 0
#
# dictionary to combine identical sequences and dictionary for pseudo gene sequences.
#
sd = defaultdict(int)
sd_pseudo = defaultdict(int)
print("Step 2: Primer check and pseudogene check!")
#
# read fastq file, check Q, check primers, combine indentical
# Check for pseudogene
#
with open(inFile,'r') as fh:  
  while True:
    id = fh.readline().strip()
    if id == "": break
    seq = fh.readline().strip()
    p   = fh.readline().strip()
    q   = fh.readline().strip()
    fail = checkQuality(q)
    if fail >= 4: 
      failQual += 1
    else:
      trim = checkPrimers(seq,start_primer,end_primer,refseq)
      if trim == "":
        failPrimers += 1
      else:
        # Pseudogene check for CD33E1 
        pgene = trim.find(pseudo[:6])
        # Pseudogene check for BE
        #pgene = trim.find("ACGACAAGAACTCCCC")
        if pgene == -1:
            sd[trim] += 1
            goodCount += 1        
        else:
            sd_pseudo[trim[54:181]] += 1
            pseudoCount += 1
    allCount += 1
#
#
#
print("total count  = {}".format(allCount))
print("fail quality = {}".format(failQual))
print("fail primers = {}".format(failPrimers))
print("total good   = {}".format(goodCount))
print("pseudo genes = {}".format(pseudoCount))
print(float(goodCount)/float(allCount))

#
# final output
#
#
print("Wild type analysis!")
needle_out = needle(sd,outFile,refFile,samFile)

print("Step 4: Identifyting substitutions, insertions and deletions!")
subs_out = seq_group(samFile,needle_out,refseq)

print("Step 5: Final step!")
out = alignment(csvFile,mutFile,subs_out,goodCount,refseq)

# Results for pseudogene
print("Pseudogene analysis!")
pseudo_needle_out = needle(sd_pseudo,pseudoFile,refFile2,pSamFile)

print("Step 4: Identifyting substitutions, insertions and deletions in pseudogenes!")
pseudo_subs_out = seq_group(pSamFile,pseudo_needle_out,pseudo)

print("Step 5: Final step for pseudogenes!")
out_pseudo = alignment(pcsvFile,pmutFile,pseudo_subs_out,pseudoCount,pseudo)




