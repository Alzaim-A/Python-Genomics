"""This script performs UPGMA hierarchical clustering on a set of species using a given distance matrix.
It iteratively merges the closest clusters using a weighted average approach. It then generates a dendrogram 
using average linkage to visualize the hierarchical clustering.
"""
import random
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

distanceMatrix_tree = [
    [0, 12, 12, 13, 15, 15],
    [12, 0, 2, 6, 8, 8],
    [12, 2, 0, 6, 9, 9],
    [13, 6, 6, 0, 8, 8],
    [15, 8, 9, 8, 0, 4],
    [15, 8, 9, 8, 4, 0]
]
speciesList_tree = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]

#copy the original distance matrix and species list for use in the UPGMA function
#save the original distance matrix and species list to create the dendrogram
speciesList = speciesList_tree.copy()          
distanceMatrix = [row.copy() for row in distanceMatrix_tree]  

def findSmallest(dM):
    """Find the indices of the smallest off-diagonal element in the distance matrix.
    param: dM - 2D list representing the symmetric distance matrix
    return: (i, j) tuple where i and j are the row and column indices of the smallest element
    """
    smallest = float('inf')
    candidates = []
    n = len(dM)
    #loop through upper triangle to find the smallest distance
    for i in range(n):
        #avoid diagonal 0's
        for j in range(i + 1, n):
            if dM[i][j] < smallest:
                smallest = dM[i][j]
                candidates = [(i, j)]
            elif dM[i][j] == smallest:
                candidates.append((i, j))
    #randomly select one of the smallest distances if there are multiple
    row, col = random.choice(candidates)
    return row, col

def updateMatrix(dM, row, col, count):
    """Merge two clusters in the distance matrix using a weighted average.
    param: dM - 2D list representing the current distance matrix.
    param: row - integer index for the first cluster to merge.
    param: col - integer index for the second cluster to merge.
    param: count - list of integers representing the number of original species in each cluster.
    return: (newMat, weighted_count) newMat is the updated matrix and weighted_count is the number of species in each cluster
    """
    n = len(dM)
    newMat = []
    #initialize and copy list to hold the species count for wieghted average
    weighted_count = count[:]

    #iterate over matrix skipping rows and col being merged
    for i in range(n):
        if i == col:
            continue
        #create matrix row by row 
        new_row = []
        for j in range(n):
            if j == col:
                continue 
            if i == row and j == row:
                #create diagonal 0 
                new_row.append(0)
            elif i == row:
                #calc weighted average distance merged cluster
                weighted_distance = (count[row] * dM[row][j] + count[col] * dM[col][j]) / (count[row] + count[col])
                new_row.append(weighted_distance)
            elif j == row:
                weighted_distance = (count[row] * dM[i][row] + count[col] * dM[i][col]) / (count[row] + count[col])
                new_row.append(weighted_distance)
            else:
                new_row.append(dM[i][j])
        newMat.append(new_row)

    weighted_count[row] = count[row] + count[col]
    del weighted_count[col]

    return newMat, weighted_count

def updateSpecies(sp, r, c, dist):
    """Merge species (or clusters) at indices r and c into a string with distances and branch lengths.
    param: sp - list of species or cluster labels.
    param: r - integer index for the first cluster.
    param: c - integer index for the second cluster.
    param: dist - total distance between the two clusters in matrix
    return: sp - updated species list with the merged cluster.
    """
    #unsure if total distance was desired or the branch length of the new species to merged cluster, so both are included
    #TD = total distance, BL = branch length
    sp[r] = "(" + sp[r] + "," + sp[c] + ": TD:" + str(dist) + ": BL:" + str(dist / 2)  + ")"
    del sp[c]
    return sp

def UPGMA(dM, sp, outFile):
    """Perform UPGMA clustering on distance matrix and species list.
    param: dM - 2D list representing the initial distance matrix.
    param: sp - list of species names corresponding to the matrix rows/columns.
    param: outFile - txt file output
    return: final clustering tree.
    """
    #initialize a list to count the number of species in each cluster
    count = [1] * len(sp)

    #loop until a final 2x2 matrix is produced
    while len(dM) > 2:
        #call findSmallest and index matrix to find distance
        r, c = findSmallest(dM)
        distance = dM[r][c]

        #update species list
        sp = updateSpecies(sp, r, c, distance)

        #update the distance matrix and species count
        dM, count = updateMatrix(dM, r, c, count)

        #print output and write to txt
        out_line = "Distance matrix:\n"
        print(out_line)
        outFile.write(out_line)
        for row_data in dM:
            line = str(row_data) + "\n"
            print(row_data)
            outFile.write(line)
        out_line = "Species list: " + str(sp) + "\n"
        print(out_line)
        outFile.write(out_line)
        out_line = "---------------\n"
        print(out_line)
        outFile.write(out_line)

    return sp


with open("aalzaim_UPGMA_output.txt", "w") as outFile:
    UPGMA(distanceMatrix, speciesList, outFile)
    final_line = "Final 2x2 Matrix and Species List with Distance and Branch Length are above:\n"
    print(final_line)
    outFile.write(final_line)

#prepare the distance matrix for the dendrogram
sqr_form_dist = squareform(distanceMatrix_tree) * 0.5

#linkage matrix for the dendrogram
Z = linkage(sqr_form_dist, method='average')

plt.figure()
dendrogram(Z, labels=speciesList_tree)
plt.title("UPGMA Tree")
plt.ylabel("Branch Length")
plt.show()
