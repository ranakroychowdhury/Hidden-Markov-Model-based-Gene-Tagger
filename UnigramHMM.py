# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:45:08 2019

@author: Ranak Roy Chowdhury
"""
import operator

def emission(O, I):
    total = sum(O.values())
    O = {k: v / total for k, v in O.items()}
    total = sum(I.values())
    I = {k: v / total for k, v in I.items()}
    return O, I

def readModelFile():
    filename = "newgene.counts"
    O = dict()
    I = dict()
    i = 0
    with open(filename) as f:
        for line in f:
            if i < 7816:
                line = line.split()
                if line[2] == 'O':
                    O[line[3]] = int(line[0])
                else:
                    I[line[3]] = int(line[0])
            i += 1
    O, I = emission(O, I)
    return O, I
    
def readDevFile():
    filename = "gene.dev"
    words = []
    with open(filename) as f:
        for line in f:
            line = line.split()
            if line:
                words.append(line[0])
            else:
                words.append(' ')
    return words

def printDictionary(dic):
    sorted(dic.items(), key=operator.itemgetter(1))
    for k,v in sorted(dic.items(), key=operator.itemgetter(1))[:10]:
        print (k,v)
        
def replace(word):
    return '_RARE_'
    
def unigramHMM(words, O, I):
    tag_sequence = []
    for word in words:
        if word == ' ':
            tag_sequence.append(' ')
            continue
        else:
            if word not in O and I:
                word = replace(word)
            
            if word in O:
                val_O = O[word]
            else:
                val_O = 0
            if word in I:
                val_I = I[word]
            else:
                val_I = 0
            
            if val_O > val_I:
                tag_sequence.append('O')
            else:
                tag_sequence.append('I-GENE')
    return tag_sequence
            
def writeResult(words, tag_sequence):
    file = open("gene_dev.p1.out", "a")
    i = 0
    for word in words:
        if word == ' ':
            file.write(' ')
        else:
            file.write(word)
        file.write(' ')
        file.write(tag_sequence[i])
        file.write('\n')
        i += 1
    file.close()

def matchKeys(keys, tag_sequence):
    count = 0
    i = 0
    mismatch = []
    for key in keys:
        if key != tag_sequence[i]:
            count += 1
            mismatch.append(i)
        i += 1
    return count, mismatch
 
def readKeys():
    filename = "gene.key"
    keys = []
    with open(filename) as f:
        for line in f:
            line = line.split()
            if line:
                keys.append(line[1])
            else:
                keys.append(' ')
    return keys
    
if __name__ == "__main__":
    O, I = readModelFile()
    words = readDevFile()
    #print(words[20:30])
    tag_sequence = unigramHMM(words, O, I)
    #print(tag_sequence[0:50])
    #writeResult(words, tag_sequence)
    keys = readKeys()
    count, mismatch = matchKeys(keys, tag_sequence)    
    print(count)
    