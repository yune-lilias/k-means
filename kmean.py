from Bio import SeqIO
import numpy as np
from sklearn.manifold import MDS
from matplotlib import pyplot as plt
from sklearn.metrics.pairwise import euclidean_distances
from matplotlib.pyplot import cm
import random
import warnings

warnings.filterwarnings("ignore")

def kmeans(X,D,k,thread):
    centres = list(range(k))
    distances = [[D[i,j] for i in centres] for j in range(len(X))]
    clusters = np.array([d0.index(min(d0)) for d0 in distances])       
    cl = [[ct for ct in range(len(clusters)) if clusters[ct]==i] for i in centres]
    X_sorted = np.array([np.array([X[xx] for xx in cl[i]]) for i in range(k)])
    centres = [np.array([np.mean(X_sorted[i][:,0]),np.mean(X_sorted[i][:,1])]) for i in range(k)]
    #print(centres)
    step = 0
    while True:
       distances = euclidean_distances(X,centres)
       clusters =  np.array([np.argmin(d0) for d0 in distances])  
       cl = [[ct for ct in range(len(clusters)) if clusters[ct]==i] for i in range(k)]     
       X_sorted = np.array([np.array([X[xx] for xx in cl[i]]) for i in range(k)])
       centres2 = []
       for X_s in X_sorted:
           if(len(X_s)>0):
               centres2.append(np.array([np.mean(X_s[:,0]),np.mean(X_s[:,1])]))
       for i in range(k):
            if(i not in clusters):
               centres2.append(centres[i])
       #print(centres2)
       #print("\n")
       step = step + 1
       for ci in range(len(cl)):
             if(len(cl[ci])<=thread*len(X)):
                  while True:
                   tmp1 = X[random.randint(0,len(X)-1)]
                   if(np.linalg.norm(tmp1-centres2[ci])!=0):
                      centres2[ci] = tmp1
                      break
       if(step >=150):
          #print("clusters: \n")
          #print(cl)
          #print("clusters size: \n")
          #print([len(ct) for ct in cl])
          break
       elif(sum([euclidean_distances(centres2,centres)[i][i] for i in range(len(centres))]) < 1):
          #print("clusters: \n")
          #print(cl)
          #print("clusters size: \n")
          #print([len(ct) for ct in cl])
          break
       else:
          centres = centres2
    return cl

fasta_sequences = SeqIO.parse(open('C:\\Users\\miond\\OneDrive\\Desktop\\py\\HW2.fas'),'fasta')

sequences = []
temp = []
for fasta in fasta_sequences:
      temp = []
      name, sequence = fasta.id, str(fasta.seq)
      temp.extend(sequence)
      sequences.append(temp)
seq = np.array(sequences)
d_matrix = (seq[:,None,:]!=seq).sum(2)
embedding = MDS(n_components=2,dissimilarity="precomputed")
X = embedding.fit_transform(d_matrix)

clusters = kmeans(X,d_matrix,k,0.1)
X_clustered = [np.array([X[i] for i in clusters0]) for clusters0 in clusters]
color= cm.rainbow(np.linspace(0,1,k))
for i,c in zip(range(k),color):
   plt.scatter(X_clustered[i][:,0], X_clustered[i][:,1],color=c)
plt.show()