# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 21:18:51 2019

@author: Ranak Roy Chowdhury
"""
import operator

def readFile():
    filename = "gene.train"
    words = []
    tags = []
    with open(filename) as f:
        for line in f:
            line = line.split()
            if line:
                words.append(line[0])
                tags.append(line[1])
            else:
                words.append(' ')
                tags.append(' ')
    return words, tags

def wordCount(word, train):
    if word in train:
        train[word] += 1.0
    else:
        train[word] = 1.0
    
def countWords(words):
    train = dict()
    for word in words:
        wordCount(word, train)
    return train
    
def printDictionary(dic):
    sorted(dic.items(), key=operator.itemgetter(1))
    for k,v in sorted(dic.items(), key=operator.itemgetter(1))[:10]:
        print (k,v)

def replace(word):
    return '_RARE_'
 
def buildNewTrain(words, tags, train, threshold_frequency):
    file = open("newgene.train", "a")
    i = 0
    count = 0
    for word in words:
        if (train[word] >= threshold_frequency):
            file.write(word)
        else:
            file.write(replace(word))
            if(word.isupper()):
                count += 1
        file.write(' ')
        file.write(tags[i])
        file.write('\n')
        i += 1
    print(count)
    file.close()
    
def replaceRareWords(words, tags, threshold_frequency):
    train = countWords(words) #count frequency of current words
    '''
    #printDictionary(train)
    upper = 0
    digit = 0
    count = 0
    for key in train:
        if(train[key] < 5):
            if key[0].isupper():
                upper += 1
            if key[-1].isdigit():
                digit += 1
            count += 1
    '''
    buildNewTrain(words, tags, train, threshold_frequency) #write new training data to file

''' 
def readCounts():
    filename = "newgene.count"
    words = []
    tags = []
    with open(filename) as f:
        for line in f:
            line = line.split()
            if line:
                words.append(line[0])
                tags.append(line[1])
            else:
                words.append(' ')
                tags.append(' ')
    return words, tags
'''
    
if __name__ == "__main__":
    words, tags = readFile()
    threshold_frequency = 5
    replaceRareWords(words, tags, threshold_frequency)    
    #readCounts()
