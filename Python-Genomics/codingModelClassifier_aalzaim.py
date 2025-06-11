import math

modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']  

def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")   
    noncodingMatrix = getProbs("./noncodingModel.tab") 
    

    id2ancestorSeq = getSeq("./Ancestor.fa") 
    id2spaciiSeq = getSeq("./Spacii.fa") 
    allID = list(id2ancestorSeq.keys()) 

    with open("aalzaim_Model_output.txt", "w") as outFile:

        for ID in allID:
            cScore = 0 
            nScore = 0 
            
            #loop through codons in input files
            for i in range(0, len(id2ancestorSeq[ID]), 3):
                a_codon = id2ancestorSeq[ID][i:i+3]
                s_codon = id2spaciiSeq[ID][i:i+3]

                #find the indicies for each codon in modelCodons
                a_index = modelCodons.index(a_codon)
                s_index = modelCodons.index(s_codon)

                #use indicies to find probs and calculate the log-odds
                cScore += math.log(codingMatrix[a_index][s_index])
                nScore += math.log(noncodingMatrix[a_index][s_index])

            if cScore > nScore:
                print (ID + "is coding" +str(cScore)+" " +str(nScore))
                outFile.write(f"{ID} is coding {cScore} {nScore}\n")
            else:
                print (ID + "is NOT coding"+ str(cScore)+" "+ str(nScore))
                outFile.write(f"{ID} is NOT coding {cScore} {nScore}\n")


def getProbs(f1):
    f = open(f1)
    pMatrix = []
    for line in f:
        tmp = line.rstrip().split("\t")
        tmp = [float(i) for i in tmp]
        pMatrix.append(tmp)
    return pMatrix

def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = (line[1:].split("|")[0])
            id2seq[currkey] = ""
        else:
            id2seq[currkey] = id2seq[currkey] + line.rstrip()
    f.close()
    return id2seq
scoreModels()