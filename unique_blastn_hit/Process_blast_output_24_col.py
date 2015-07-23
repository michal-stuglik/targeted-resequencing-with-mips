'''
This script takes a file of blast results in 24 col NCBI format
Creates a dictionary; each dictionary entry contains all lines reporting hits to a given query
Within this first order dictionary a second order dictionary is created, 
containing all lines which are hsps for a given subject 
The second order dictionaries contain one (in case of a single hsp) or more lists, which
accommodate information from individual blast output_columns
Then, for each query, the number of subjects with hits is calculated and three cases are outputed into three files
1. Queries with a single subject - all hsps writen
2. Queries with multiple subjects - for each subject a "total" bitscore is calculated and:
- if the subject with the highest bitscore is much better than the next (customizable through NextBestThreshold), all hsps
  for this subject are outputted
- if the best is not much better than the best, then for all subjects with bitscore constituting at least NextBestThreshold 
  fraction of the best are outputted
USAGE: enter the name of the file containing your 24 column blast output in InputFilename below
output consists of three files
Single_hit.txt  
Multiple_several_good.txt
Multiple_one_the_best
'''

import sys

#InputFilename="output_blast"
InputFilename=sys.argv[1]


def readBlastOutputDicionary(InputFilename):
    InFile=open(InputFilename, 'r')
    QuerColl={}
    for Line in InFile:
        Line=Line.strip()
        LIL=Line.split()
        Quer=LIL[0]
        Subj=LIL[1]
        if Quer in QuerColl:				#check whether there is the main dictionary entry for the query
            SubjDict=QuerColl[Quer]			#if yes, put this entry to a temporary dictionary
            if Subj in SubjDict:			#check whether there is an entry in the temporary dictionary
                                            #(this dictionary equals an entry in the higher level dictionary) for the subject
                SubjList=SubjDict[Subj]		# if yes, put this entry to a "temporary" list
                SubjList.append(LIL)	#add the required part of the line to the list (this is another hsp for the subject contig
                SubjDict[Subj]=SubjList		#update temporary dictionary
            else:							#if there is no entry in temporary dictionary, it means that the contig is new for our query
                SubjList=[]					#create empty list
                SubjList.append(LIL)	#add the required part of the line to the list (this is the first hsp for the subject contig
                SubjDict[Subj]=SubjList		#update temporary dictionary - create its first item
            QuerColl[Quer]=SubjDict			#update the main dictionary entry for the query
        else:						#if there is no main dictionary entry for the query
            SubjDict={}				#create empty temporary dictionary
            SubjList=[]				#create empty list
            SubjList.append(LIL)		#add the required part of the line to the list (this is the first hsp for the first subject contig)
            SubjDict[Subj]=SubjList			#update temporary dictionary - create its first item
            QuerColl[Quer]=SubjDict 		#create the main dictionary item for the query, containing a single contig with single hsp
    InFile.close()
    return QuerColl
    
def writeHspToFile(HspToWrite,FileToWrite):
    OutFile=open(FileToWrite, 'a')
    for item in HspToWrite[0:12]:				#writes only the "short" table
        OutFile.write("%s\t" % item)
    OutFile.write("\n")
    OutFile.close()

def extractHspInfo(Hsp, NumberOfHsp):
    ListToAdd=[]
    for i in xrange(NumberOfHsp):
        HspL=float(Hsp[i][3])				#Hsp length
        HspBS=float(Hsp[i][11])				#extract Hsp bitscore
        MeanBS=HspBS/HspL					#calculate mean bitscore
        start=Hsp[i][6]
        end=Hsp[i][7]
        if start < end:
            HspStart=int(Hsp[i][6])
            HspEnd=int(Hsp[i][7])
        else:
            HspStart=int(Hsp[i][7])
            HspEnd=int(Hsp[i][6])
        ListToAdd.append([HspStart, HspEnd, MeanBS])
    return ListToAdd

def calculateSubjectBitscore(LengthOfQuery,NumberOfHsp,ListOfParams):
    SubjectBitscore=0
    for base in range(1, LengthOfQuery+1):
        BaseBitscore=0
        HspCounter=0
        for i in range(NumberOfHsp):
            HspStart=ListOfParams[i][0]
            HspEnd=ListOfParams[i][1]
            MeanBS=ListOfParams[i][2]
            if (HspStart<=base and HspEnd>=base):
                BaseBitscore = BaseBitscore + MeanBS
                HspCounter=HspCounter + 1
            else:
                    pass
        if HspCounter==0:
            pass
        else:
            SubjectBitscore=SubjectBitscore + float(BaseBitscore/HspCounter)
    return SubjectBitscore

def keywithmaxval(d):
#a) create a list of the dict's keys and values; 
#b) return the key with the max value"""  
    v=list(d.values())
    k=list(d.keys())
    return k[v.index(max(v))]



#BODY
BlastRes=readBlastOutputDicionary(InputFilename)
NextBestThreshold=0.8
for Query in BlastRes:									#process each query
    NSubj=len(BlastRes[Query])
    if NSubj==1:
        for Subject in BlastRes[Query]:					#this could be probably put into a short general function
            for Hsp in BlastRes[Query][Subject]:
                writeHspToFile(Hsp,'Single_hit.txt')
    else:
        BitscoreDict={}									#Keeps bitscores for subjecs
        for Subject in BlastRes[Query]:					#process each subject
            NHsp=len(BlastRes[Query][Subject])			#how many hsp?
            Hsps=BlastRes[Query][Subject]
            if NHsp==1:
                BitscoreDict[Subject]=float(Hsps[0][11]) #if single hsp gets its bitscore
            else:
                ListOfHspParams = extractHspInfo(Hsps, NHsp)	# extracts Hsp's per base bitscore, start and end
                QueryL=int(Hsps[0][22])							# calculate query length using data for first Hsp (doesn't matter which one is used)
                SubjectBitscore = calculateSubjectBitscore(QueryL, NHsp, ListOfHspParams) #calclulates SubjectBitscore in case of multiple Hsps
                BitscoreDict[Subject]=SubjectBitscore									#puts subject's bitscore into a dictionary
        BitscoreList=[]
        for key, value in BitscoreDict.iteritems():										#converts dictionary into list
            temp=[key,float(value)]
            BitscoreList.append(temp)
        SortedBL=sorted(BitscoreList, key=lambda Bitscore: Bitscore[1], reverse=True)	  #sorts the list of subjects descending accoring to bitscore
        if SortedBL[1][1] < SortedBL[0][1]*NextBestThreshold:			# if best is best enough then retain only the best subject
            for Hsp in BlastRes[Query][SortedBL[0][0]]:
                writeHspToFile(Hsp,'Multiple_one_the_best.txt')
        else:																				#if not best enough output all with bitscores at least NextBestThreshold*best
            i=0
            while ((i+1)<=len(SortedBL) and (SortedBL[i][1] >= (SortedBL[0][1]*NextBestThreshold))):
                for Hsp in BlastRes[Query][SortedBL[i][0]]:
                    writeHspToFile(Hsp,'Multiple_several_good.txt')
                i = i + 1
