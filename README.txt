***[CSE 256 FA 19: Programming assignment 3: Sequence Tagging]***

***[ Files ]***

There are five python files in this folder and the README file:

-(baseline.py): This file reads the training data from gene.train. Replaces words that appear less than 5 times by informative classes and builds the new training data newgene.train.

- (count_freqs.py): This file builds the count file from the training data. The count file contains the number of times each word appear with the tag 'O' and 'I-GENE'. It also contains the frequency of all posiible trigrams, bigrams and unigrams. 

You should run: python count_freqs.py newgene.train > newgene.counts

The above command produces the count file from the training data.

- (UnigramHMM.py): This file contains the Baseline Tagger. Reads the newgene.counts and the gene.dev file to produce the gene_dev.p1.out file that contains the words from the dev file and the predicted tag associated with each word.

- (TrigramHMM.py): This file contains the Trigram HMM Tagger. Reads the newgene.counts and the gene.dev file to produce the gene_dev.p2.out file that contains the words from the dev file and the predicted tag associated with each word.

- (eval_gene_tagger.py): This file compares the gene.key file with the file produced by the tagger. It generates the precision, recall and F1 score.

You should run: python eval_gene_tagger.py gene.key gene_dev.p1.out

The above command generates the performance measures for the UnigramHMM Tagger.

You should run: python eval_gene_tagger.py gene.key gene_dev.p2.out

The above command generates the performance measures for the TrigramHMM Tagger.


