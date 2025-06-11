import math
import matplotlib.pyplot as plt  

baseIDx = {"A":0, "C":1, "G":2, "T":3}

def main():
    spaciiFA = "MSpacii.fa"
    pathogenFA = "pathogen.fa"
    spaciiFA_T = "MSpacii_training.fa"
    pathogenFA_T = "pathogen_training.fa"
    spaciiID2seq = getSeq(spaciiFA)
    pathogenID2seq = getSeq(pathogenFA)
    
    spaciiTrainModel = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    pathTrainModel = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    
    spaciiTrainModel = trainModel(spaciiTrainModel, spaciiFA_T)
    pathTrainModel = trainModel(pathTrainModel, pathogenFA_T)
    
    markovScoresSpacii = []
    markovScoresPath = []

    for ID in spaciiID2seq.keys():
         markovScoresSpacii.append(getLogLike(spaciiTrainModel, pathTrainModel, spaciiID2seq[ID]))
    for ID in pathogenID2seq.keys():
         markovScoresPath.append(getLogLike(spaciiTrainModel, pathTrainModel, pathogenID2seq[ID]))
    
    ####----------------------output------------------------- 
    #commented out per notes   
    plt.hist([markovScoresPath, markovScoresSpacii], bins=20, label=['Pathogen','Spacii'], rwidth=1, density=False)
    plt.title("Distribution of Pathogen and Spacii Scores")
    plt.legend(loc='upper right')
    plt.xlabel("Log-Likelihood")
    plt.ylabel("Frequency")

    plt.show()   
    scoresOutputText(markovScoresSpacii, markovScoresPath)
    ####----------------------output-------------------------

def scoresOutputText(markovScoresSpacii, markovScoresPath):
    f = open("results.tab", "w")
    f.write("SpaciiScores\tpathogenScores\n")
    for i in range(len(markovScoresSpacii)):
        f.write(str(markovScoresSpacii[i]) + "\t" + str(markovScoresPath[i]) + "\n")
    f.close()
    
def getLogLike(model1, model2, seq):
    
    Pmod1 = 1
    Pmod2 = 1
    
    #iterate through each base and the previousl
    for i in range(1, len(seq)):
        prev_base = seq[i - 1]
        curr_base = seq[i]
        idx_prev = baseIDx[prev_base]
        idx_curr = baseIDx[curr_base]
        #get the probability of the dinucleotide pair from each model
        p1 = model1[idx_prev][idx_curr]
        p2 = model2[idx_prev][idx_curr]
        
        #log
        Pmod1 += math.log(p1)
        Pmod2 += math.log(p2)
    
    #log-likelihood ratio
    score = Pmod1 - Pmod2
    return score
    
def trainModel(model, data):
    #get sequences from the fasta file
    id2seq = getSeq(data)
    
    #interate through each base in sequence and the previous base to get counts
    for seq in id2seq.values():
        for i in range(1, len(seq)):
            prev_base = seq[i - 1]
            curr_base = seq[i]
            model[baseIDx[prev_base]][baseIDx[curr_base]] += 1
            
    #convert raw counts to probabilities
    for i in range(4):
        row_sum = sum(model[i])
        if row_sum > 0:
            model[i] = [count / row_sum for count in model[i]]

    print(data)
    print(model)
    return model  

def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = line.rstrip()[1:]
            id2seq[currkey] = ""
        else:            
            id2seq[currkey] = id2seq[currkey] + line.rstrip()
    return id2seq

main()
