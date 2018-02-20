# Bayesian Gene Regulatory Network

**Author**:		   Danielle Novick  
**Course**:		   CISC889 - Modeling and Simulation in Bioinformatics  
**Objective**:  Reconstruct a gene regulatory network from binary gene expression data (extremely small mock data) by implementing Bayesian networks and make inferences about gene expression  
**Location**:		https://github.com/juniper-lake/Bayesian_Gene_Regulatory_Network.git
            
         
### Running Program

* To run program with default settings:  
`$ python network_generator.py`

* To show help file:  
`$ python network_generator.py --help`

* Example of modified parameters using long names:
`$ python network_generator.py --testGene 3 --maxParents 4 --maxAttempts 100`

* Example of modified parameters using short names:
`$ python network_generator.py -g 3 -p 4 -a 100`



### Help File 

``` 
$ python network_generator.py --help
usage: network_generator.py [-h] [--testGene TESTGENE]
                               [--maxParents MAXPARENTS]
                               [--maxAttempts MAXATTEMPTS]
                               [--trainingData TRAININGDATA]
                               [--testData TESTDATA]

optional arguments:
  -h, --help            show this help message and exit
  --testGene TESTGENE, -g TESTGENE
                        The gene that will be predicted in test data to report
                        prediction accuracy with first gene number starting at
                        0 so options using this data are (0,1,2,3,4,5),
                        default is gene 4
  --maxParents MAXPARENTS, -p MAXPARENTS
                        The maximum number of parent nodes that each node can
                        have, default is 3 parents
  --maxAttempts MAXATTEMPTS, -a MAXATTEMPTS
                        The maximum number of graphs that will be proposed
                        without improving upon the existing best graph before
                        program ends, default is 200
  --trainingData TRAININGDATA, -train TRAININGDATA
                        The source file of the training data, must be in same
                        format as was provided for HW, default is train.txt
  --testData TESTDATA, -test TESTDATA
                        The source file of the test data, must be in same
                        format as was provided for HW, default is test.txt
```

### Example Output

```
Generating initial graph...
Initial score:  -340.906166444
Initial graph:  [(0, 5), (2, 4), (3, 4), (5, 4)] 

Updated score:  -336.674907206 		Attempts:  3
Updated graph:  [(0, 5), (2, 3), (2, 4), (3, 4), (5, 4)] 

Updated score:  -336.674907206 		Attempts:  4
Updated graph:  [(2, 3), (2, 4), (3, 4), (5, 0), (5, 4)] 

Updated score:  -336.4548583 		Attempts:  1
Updated graph:  [(2, 0), (2, 3), (2, 4), (3, 4), (5, 0), (5, 4)] 

Updated score:  -336.453393037 		Attempts:  1
Updated graph:  [(2, 0), (2, 3), (2, 4), (3, 4), (5, 0), (5, 2), (5, 4)] 

Updated score:  -336.453393037 		Attempts:  3
Updated graph:  [(0, 2), (2, 3), (2, 4), (3, 4), (5, 0), (5, 2), (5, 4)] 

Updated score:  -333.804596607 		Attempts:  4
Updated graph:  [(0, 2), (2, 3), (2, 4), (3, 4), (5, 0), (5, 1), (5, 2), (5, 4)] 

Updated score:  -332.840879615 		Attempts:  1
Updated graph:  [(0, 2), (1, 2), (2, 3), (2, 4), (3, 4), (5, 0), (5, 1), (5, 2), (5, 4)] 

Updated score:  -332.8383832 		Attempts:  1
Updated graph:  [(0, 2), (1, 0), (1, 2), (2, 3), (2, 4), (3, 4), (5, 0), (5, 1), (5, 2), (5, 4)] 

Updated score:  -331.680296649 		Attempts:  10
Updated graph:  [(0, 2), (1, 0), (1, 2), (2, 3), (2, 4), (3, 4), (5, 0), (5, 1), (5, 2), (5, 3), (5, 4)] 

Updated score:  -331.680296649 		Attempts:  10
Updated graph:  [(0, 2), (1, 0), (1, 2), (1, 5), (2, 3), (2, 4), (3, 4), (5, 0), (5, 2), (5, 3), (5, 4)] 

Updated score:  -330.719991839 		Attempts:  18
Updated graph:  [(0, 2), (1, 0), (1, 2), (1, 3), (1, 5), (2, 3), (2, 4), (3, 4), (5, 0), (5, 2), (5, 3), (5, 4)] 

Final score:  -330.719991839 		Attempts:  200
Final graph:  [(0, 2), (1, 0), (1, 2), (1, 3), (1, 5), (2, 3), (2, 4), (3, 4), (5, 0), (5, 2), (5, 3), (5, 4)] 

Prediction accuracy for gene 4:  0.2 

Conditional probability tables of final graph: 

Gene:  0 
    1  5       Pr0       Pr1
0  0  0  0.296875  0.703125
1  0  1  0.272727  0.727273
2  1  0  0.294118  0.705882
3  1  1  0.264151  0.735849 

Gene:  1 
     Pr0   Pr1
0  0.65  0.35 

Gene:  2 
    0  1  5       Pr0       Pr1
0  0  0  0  0.684211  0.315789
1  0  0  1  0.555556  0.444444
2  0  1  0  0.200000  0.800000
3  0  1  1  0.571429  0.428571
4  1  0  0  0.466667  0.533333
5  1  0  1  0.500000  0.500000
6  1  1  0  0.583333  0.416667
7  1  1  1  0.487179  0.512821 

Gene:  3 
    1  2  5       Pr0       Pr1
0  0  0  0  0.617647  0.382353
1  0  0  1  0.617647  0.382353
2  0  1  0  0.433333  0.566667
3  0  1  1  0.281250  0.718750
4  1  0  0  0.875000  0.125000
5  1  0  1  0.703704  0.296296
6  1  1  0  0.666667  0.333333
7  1  1  1  0.230769  0.769231 

Gene:  4 
    2  3  5       Pr0       Pr1
0  0  0  0  0.285714  0.714286
1  0  0  1  0.550000  0.450000
2  0  1  0  0.357143  0.642857
3  0  1  1  0.285714  0.714286
4  1  0  0  0.421053  0.578947
5  1  0  1  0.466667  0.533333
6  1  1  0  0.150000  0.850000
7  1  1  1  0.325581  0.674419 

Gene:  5 
    1       Pr0       Pr1
0  0  0.492308  0.507692
1  1  0.242857  0.757143 
```
