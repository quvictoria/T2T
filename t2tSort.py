# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 22:48:28 2020

@author: Victoria Qu
"""

#import all libraries needed
import  urllib.request
import gzip
import numpy as np 
import pandas as pd
import csv

#main method to call other methods 
def main():
    #(summarizeData('myOutFile.txt'))
    print(determine(summarizeData('myOutFile.txt')))
    data = determine(summarizeData('myOutFile.txt'))
    dataFrame = pd.DataFrame(data, columns = ['Convergent Gene Pair',
                                              'Base Pairs Apart',
                                              'Coordinates of Gene A (+ strand)',
                                              'Coordinates of Gene B (- strand)', 
                                              'Location'])
    dataFrame.to_csv('t2tORF.csv')
    
#retrieves data from the SGD website 
"""
What it does: retrieves a .gz file from the SGD archives by sending a url library request, 
unzips the contents, converts to bytes and retrieves relevant data to be printed to a .txt file which 
will be used later. 
What it returns: doesnt return anything, just creates 'myOutFile.txt'
"""
def retrieveData():
    url='http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz'
    r =  urllib.request.urlretrieve(url, 'my.gz')
    outF = open("myOutFile.txt", "wb")
    with gzip.open('my.gz', 'rb') as f:
        for line in f:
            if(line[0:3]==b'chr'):
                txt=line.split(b'\t')
                if (len(txt)>7):
                    if(txt[2]==b'ARS' or txt[2]==b'telomere' or txt[2]==b'centromere' or (b'Ontology_term=' in txt[8])):
                       flg=True               
                       cr=txt[0]
                       tp=txt[2]
                       if(tp==b'gene'):
                           tp=b'ORF'
                       loc1=txt[3]
                       loc2=txt[4]
                       stra=txt[6]
                       note=txt[8]
                       nts=note.split(b';')
                       for nt in nts:
                           ntt=nt.split(b'=')
                           if(len(ntt)==2):
                               if(ntt[0]==b'ID'):
                                   id=ntt[1]
                               if(ntt[0]==b'Parent'):
                                    flg=False
                               if(ntt[0]==b'Name'):
                                   name=ntt[1]
                                
                               if(ntt[0]==b'gene'):
                                   name=ntt[1]
                         
                       if (flg):
                          outF.write( cr+b'\t'+name +b'\t'+ id+b'\t'+ tp+b'\t'+ loc1+b'-'+loc2+b'\t'+stra+b'\n')
    f.close()
    outF.close()

#summarizes data and gets rid of nonrelevant entries 
"""
What it does: reads in file that we created from retrieval method and stores only ORFs it in a Numpy array 
What it returns: Returns a Numpy array of only ORF genes 
"""
def summarizeData(filename):
    #read in excel file into numpy array
    data = np.genfromtxt(filename, dtype = ['<U16', '<U16', '<U16', '<U16', '<U16', '<U16'], delimiter = '\t')
    
    #initialize array 
    arr = []
    #loop through array generated from .cvs file     
    for i in data:
        #only load data that is relevant to ORFs 
        if (i[3] == 'ORF'):
            temp = [i[0], i[1], i[2], i[3], i[4], i[5]]
            arr.append(temp)
    #return array        
    return np.array(arr)
"""
What it does: Loops through previous np array to find pairs and then checks the distance 
What it returns: Np array with the pairs, how far they are from each other, location
"""
def determine(arr):
    #loop through array 
    arr1 = []
    prev = arr[0][5]
    #print(prev)
    i= 1
    howFarBy = 100
    while (i < len(arr)):
        if arr[i-1][5] == '+' and arr[i][5] == '-':
            a = arr[i-1][4].split("-")
            b = arr[i][4].split("-")
            aEnd = a[1]
            bStart = b[0]
            distance = int(bStart) - int(aEnd)
            pair = arr[i-1][1] + " and " + arr[i][1]
            
            #checks distance between convergent pairs, the number can be changed by changing variable howFarBy
            if distance < howFarBy and distance>0: 
                temp = [pair, distance, arr[i-1][4], arr[i][4], arr[i][0]]
                arr1.append(temp)
            """
            #create new 
            
            pair = arr[i-1][0] + " and " + arr[i][0]
            #print(pair)
            a = arr[i-1][3].split("-")
            b = arr[i][3].split("-")
            distance = int(b[0])-int(a[1])
            if 0 <= distance <= 100:
                #array will have: pair names, distance between, coordinates of gene 1, coordinates of gene2 
                temp = [pair, distance, arr[i-1][3], arr[i][3]]
                arr1.append(temp)"""
            
        i+=1
    return np.array(arr1)
    

if __name__ == "__main__": main() 