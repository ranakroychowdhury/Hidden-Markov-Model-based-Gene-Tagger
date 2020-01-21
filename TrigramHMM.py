# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:59:26 2019

@author: Ranak Roy Chowdhury
"""
import operator

def emission(O, I):
    total = sum(O.values())
    O = {k: v / total for k, v in O.items()}
    total = sum(I.values())
    I = {k: v / total for k, v in I.items()}
    return O, I

def transitionProbability(bigram, trigram):
    dic = dict()
    for pair in bigram:
        for triplet in trigram:
            if pair[0] == triplet[0] and pair[2] == triplet[2]:
                key = triplet[0] + ' ' + triplet[2] + ' ' + triplet[4]
                dic[key] = int(triplet[5])/int(pair[3])
    return dic
    
def readModelFile():
    filename = "newgene.counts"
    O = dict()
    I = dict()
    bigram = []
    trigram = []
    i = 0
    with open(filename) as f:
        for line in f:
            line = line.split()
            if i < 7816:
                if line[2] == 'O':
                    O[line[3]] = int(line[0])
                else:
                    I[line[3]] = int(line[0])
            else:
                l = []
                if i > 7817 and i < 7827:
                    l.append(line[2])
                    l.append(' ')
                    l.append(line[3])
                    l.append(line[0])
                    bigram.append(l)
                elif i >= 7827:
                    l.append(line[2])
                    l.append(' ')
                    l.append(line[3])
                    l.append(' ')
                    l.append(line[4])
                    l.append(line[0])
                    trigram.append(l)
            i += 1
    O, I = emission(O, I)
    dic = transitionProbability(bigram, trigram)
    return O, I, dic

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
    
def backtrack(tag_sequence, tracker_O, tracker_I):
    end = len(tracker_O) - 1
    for i in range(len(tracker_O)):
        if tag_sequence[i + 1] == 0:
            tag_sequence.append(tracker_O[end])
        else:
            tag_sequence.append(tracker_I[end])
        end = end - 1
    return tag_sequence

def process(tag_sequence):
    tag_sequence = tag_sequence[::-1]
    new_tag_sequence = []
    for i in range(len(tag_sequence)):
        if tag_sequence[i] == 0:
            new_tag_sequence.append('O')
        else:
            new_tag_sequence.append('I-GENE')
    return new_tag_sequence
    
def viterbi(sentence, emission_O, emission_I, transition):
    p_O = []
    p_I = []
    tracker_O = []
    tracker_I = []
    tag_sequence = []
    
    #1st word
    word_count = 0
    if sentence[0] in emission_O:
        p_O.append(transition['* * O'] * emission_O[sentence[0]])
    else:
        p_O.append(transition['* * O'] * 0)
    if sentence[0] in emission_I:
        p_I.append(transition['* * I-GENE'] * emission_I[sentence[0]])
    else:
        p_I.append(transition['* * I-GENE'] * 0)
        
    #2nd word
    word_count = 1
    if word_count < len(sentence) - 1:
        if sentence[1] in emission_O:
            p_O.append(max(transition['* O O'] * emission_O[sentence[1]], transition['* I-GENE O'] * emission_O[sentence[1]]))
        else:
            p_O.append(max(transition['* O O'] * 0, transition['* I-GENE O'] * 0))
        if sentence[1] in emission_I:
            p_I.append(max(transition['* O I-GENE'] * emission_I[sentence[1]], transition['* I-GENE I-GENE'] * emission_I[sentence[1]]))
        else:
            p_I.append(max(transition['* O I-GENE'] * 0, transition['* I-GENE I-GENE'] * 0))
    else:
        l_O = []
        l_I = []
        final_list = []
        final_tags = []
        l_O.append(transition['* O STOP'])
        l_O.append(transition['* I-GENE STOP'])
        l_I = l_O
        if sentence[1] in emission_O:
            l_O = [i * emission_O[sentence[1]] for i in l_O]
        else:
            l_O = [i * 0 for i in l_O]
        if sentence[1] in emission_I:
            l_I = [i * emission_I[sentence[1]] for i in l_I] 
        else:
            l_I = [i * 0 for i in l_I]
        l_O_max = max(l_O)
        l_O_max_idx = l_O.index(max(l_O))
        l_I_max = max(l_I)
        l_I_max_idx = l_I.index(max(l_I))
        final_list.append(l_O_max)
        final_list.append(l_I_max)
        final_tags.append(l_O_max_idx)
        final_tags.append(l_I_max_idx)
        tag_id = final_list.index(max(final_list))
        if tag_id == 0:
            tag_sequence.append(0)
        else:
            tag_sequence.append(1)
                
        last_tag = final_tags[tag_id]
        if last_tag == 0:
            tag_sequence.append(0)
        else:
            tag_sequence.append(1)
        tag_sequence = process(tag_sequence)
        return tag_sequence
    
    #Rest of the words
    for word in sentence[2:]:
        word_count += 1
        if word_count < len(sentence) - 1:
            l_O = []
            l_O.append(p_O[word_count - 2] * transition['O O O'])
            l_O.append(p_O[word_count - 2] * transition['O I-GENE O'])
            l_O.append(p_I[word_count - 2] * transition['I-GENE O O'])
            l_O.append(p_I[word_count - 2] * transition['I-GENE I-GENE O']) 
            if word in emission_O:
                l_O = [i * emission_O[word] for i in l_O]
            else:
                l_O = [i * 0 for i in l_O]
            p_O.append(max(l_O))
            idx = l_O.index(max(l_O))
            if idx == 0 or idx == 1:
                tracker_O.append(0)
            else:
                tracker_O.append(1)
                    
            l_I = []
            l_I.append(p_O[word_count - 2] * transition['O O I-GENE'])
            l_I.append(p_O[word_count - 2] * transition['O I-GENE I-GENE'])
            l_I.append(p_I[word_count - 2] * transition['I-GENE O I-GENE'])
            l_I.append(p_I[word_count - 2] * transition['I-GENE I-GENE I-GENE'])
            if word in emission_I:
                l_I = [i * emission_I[word] for i in l_I]
            else:
                l_I = [i * 0 for i in l_I]
            p_I.append(max(l_I))
            idx = l_I.index(max(l_I))
            if idx == 0 or idx == 1:
                tracker_I.append(0)
            else:
                tracker_I.append(1)
        else:
            final_list = []
            final_tags = []
            l_O = []
            l_I = []
            l_O.append(p_O[word_count - 2] * transition['O O STOP'])
            l_O.append(p_O[word_count - 2] * transition['O I-GENE STOP'])
            l_O.append(p_I[word_count - 2] * transition['I-GENE O STOP'])
            l_O.append(p_I[word_count - 2] * transition['I-GENE I-GENE STOP'])
            l_I = l_O
            if word in emission_O:
                l_O = [i * emission_O[word] for i in l_O]
            else:
                l_O = [i * 0 for i in l_O]
            if word in emission_I:
                l_I = [i * emission_I[word] for i in l_I]   
            else:
                l_I = [i * 0 for i in l_I]
            l_O_max = max(l_O)
            l_O_max_idx = l_O.index(max(l_O))
            l_I_max = max(l_I)
            l_I_max_idx = l_I.index(max(l_I))
            final_list.append(l_O_max)
            final_list.append(l_I_max)
            final_tags.append(l_O_max_idx)
            final_tags.append(l_I_max_idx)
            tag_id = final_list.index(max(final_list))
            
            if tag_id == 0:
                tag_sequence.append(0)
            else:
                tag_sequence.append(1)
                
            last_tag = final_tags[tag_id]
            if last_tag == 0:
                tag_sequence.append(0)
                tag_sequence.append(0)
            elif last_tag == 1:
                tag_sequence.append(1)
                tag_sequence.append(0)
            elif last_tag == 2:
                tag_sequence.append(0)
                tag_sequence.append(1)
            else:
                tag_sequence.append(1)
                tag_sequence.append(1)
    
    tag_sequence = backtrack(tag_sequence, tracker_O, tracker_I)
    new_tag_sequence = process(tag_sequence)
    return new_tag_sequence

def replace(words, O, I):
    new_words = []
    for word in words:
        if word == ' ':
            new_words.append(word)
        elif word in O or word in I:
            new_words.append(word)
        else:
            new_words.append('_RARE_')
    return new_words
        
def makeSentence(words):
    sentence = []
    indices = [i for i, x in enumerate(words) if x == ' ']
    start = 0
    for i in range(len(indices)):
        temp = words[start : indices[i]]
        sentence.append(temp)
        start = indices[i] + 1
    return sentence

def writeResult(words, tag_sequence):
    file = open("gene_dev.p2.out", "a")
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
    
if __name__ == "__main__":
    O, I, transition = readModelFile()
    words = readDevFile()
    #printDictionary(transition)
    replaced_words = replace(words, O, I)
    sentences = makeSentence(replaced_words)
    
    tag = viterbi(sentences[0], O, I, transition)
    print(tag)
    
    tag_sequence = []
    for sentence in sentences:
        tags = viterbi(sentence, O, I, transition)
        tags.append(' ')
        tag_sequence.append(tags)
    tag_sequence = sum(tag_sequence, [])
    #print(len(words))
    #print(len(tag_sequence))
    writeResult(words, tag_sequence)
    