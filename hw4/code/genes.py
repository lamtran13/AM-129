import numpy as np
import matplotlib.pyplot as plt

def codonHistogram(cov2):
    d = {}
    nCodons = len(cov2)-3
    for i in range(0,nCodons):
        c = scov2[i:i+3]
        d[c] = d.get(c,0) + 1
    for c,val in d.items():
        d[c] = val/nCodons
    return d

def plotHistogram(histo, title, filename):
    plt.figure(figsize = (12,6))
    plt.bar(range(len(histo)), histo.values())
    plt.xticks(ticks = np.arange(len(histo)), labels=list(histo), rotation = 'vertical')
    plt.xlabel('Codon')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.title(str(title))
    plt.savefig(str(filename))
    plt.show()

####################################
def baseDensity(geneStr, nWind=200):
    densityA = np.zeros(len(geneStr)-nWind)
    densityC = np.zeros(len(geneStr)-nWind)
    densityT = np.zeros(len(geneStr)-nWind)
    densityG = np.zeros(len(geneStr)-nWind)

    length = len(geneStr)-nWind
    for i in range(0,length):
            densityA[i] = geneStr[i:i + nWind].count('A')/nWind
            densityC[i] = geneStr[i:i + nWind].count('C')/nWind
            densityT[i] = geneStr[i:i + nWind].count('T')/nWind
            densityG[i] = geneStr[i:i + nWind].count('G')/nWind
            
            return densityA, densityC, densityT, densityG

def plotHistogram2(histo,title, filename2, var='dens',color='black'):
    plt.figure(figsize = (12,6))
    plt.plot(densityA,color=color,linewidth=2)
    #plt.bar(range(len(histo)), histo.values())
    plt.xticks(ticks = np.arange(len(histo)), labels=list(histo), rotation = 'vertical')
    plt.title(str(title))
    plt.xlabel('Sequence position')
    plt.ylabel('Fraction per window')
    plt.grid(True)
    plt.savefig(str(filename2))
    plt.show()

if __name__=="__main__":
    with open('seqSC2.txt','r') as geneFile:
        scov2 = geneFile.readline()

    with open('seqA.txt','r') as seqA:
        geneA = seqA.readline()

    with open('seqB.txt','r') as seqB:
        geneB = seqB.readline()
    
    histSC2 = codonHistogram(scov2)

    plotHistogram(histSC2,'Codon frequency in SARS-CoV-2','histoSC2.png')

    histA = baseDensity(geneA)
    plotHistogram2(histA, 'Sequence position', 'density.png', var='dens',color='blue')

    histB = baseDensity(geneB)
    plotHistogram2(histB, 'Sequence position', 'density.png',var='dens',color='green')
    
    # Find base-pair density
    densityA, densityT, densityC, densityG = baseDensity(scov2)


