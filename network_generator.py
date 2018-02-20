
#!/usr/bin/python3

# - - - - - H E A D E R - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
'''
Author:		Danielle Novick
Course:		CISC889: Modeling and Simulation in Bioinformatics
Assignment:	Homework 2
Date Due:	17 October 2017

Objective:  Reconstruct a gene regulatory network from gene expression data by implementing Bayesian networks
            and make inference about gene expression.

'''

import random
import numpy as np
import itertools
import pandas as pd
from collections import defaultdict
import argparse


# - - - - - U S E R   V A R I A B L E S - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - G L O B A L  D E C L A R A T I O N S  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
parser = argparse.ArgumentParser()

parser.add_argument('--testGene','-g', action="store", default=4, type=int,
                    help='The gene that will be predicted in test data to report prediction accuracy with first gene number starting at 0 so options using this data are (0,1,2,3,4,5), default is gene 4')
parser.add_argument('--maxParents','-p', action="store", default=3, type=int,
                    help='The maximum number of parent nodes that each node can have, default is 3 parents')
parser.add_argument('--maxAttempts', '-a', action="store", default=200, type=int,
                    help='The maximum number of graphs that will be proposed without improving upon the existing best graph before program ends, default is 200')
parser.add_argument('--trainingData', '-train', action="store", default="train.txt", dest="trainingData", type=str,
                    help='The source file of the training data, must be in same format as was provided for HW, default is train.txt')
parser.add_argument('--testData', '-test', action="store", default="test.txt", dest="testData", type=str,
                    help='The source file of the test data, must be in same format as was provided for HW, default is test.txt')

args = parser.parse_args()
testGene = args.testGene
maxParents = args.maxParents
maxAttempts = args.maxAttempts
trainingData = args.trainingData
testData = args.testData


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - P A R T   1 - - - - -
# Regulatory Network Construction:  Implement the learning algorithm based on maximum likelihood.

def read_data(filename):
    """
    Import data as a dataframe
    :param filename: name of file to be imported
    :return: a dataframe
    """
    data = pd.DataFrame.from_csv(filename, sep='\t', header=None)
    data = data.reset_index(drop=True)
    data.columns = list(range(len(data.columns)))
    return data


def edges_to_graph(edges):
    """
    Convert a list of edges to a dictionary
    :param edges: a list of tuples
    :return: a dictionary with parent nodes as keys and child nodes as values
    """
    graph = defaultdict(list)
    edges.sort()
    for key, value in edges:
        graph[key].append(value)
    return graph


def random_directed_graph(data, maxParents):
    """
    Build a random directed graph.
    :param data: the training data that will be used with the graph
    :param maxParents: the maximum number of parents for each node
    :return: a random directed graph in the form of a dictionary
    """
    numNodes = len(data.columns)
    edges = []
    nodes = list(range(numNodes))
    for node in nodes:
        otherNodes = list(range(0, node)) + list(range(node+1, numNodes))
        numEdges = random.choice(range(maxParents+1))
        for edge in range(numEdges):
            random.shuffle(otherNodes)
            parent = otherNodes.pop()
            edges.append((parent, node))
    edges.sort()
    return edges


def isCyclicUtil(edges, v, visited, recStack):
    """
    Utility function for isCyclic()
    Function is modified from http://www.geeksforgeeks.org/detect-cycle-in-a-graph/
    :param graph: a directed graph in the form of a dictionary
    :param parent: each key in the dictionary "graph"
    :param visited: a list that keeps track of which nodes have been checked already
    :param recStack: recursion stack to build paths to all potential neighbors
    :return: boolean
    """
    graph = edges_to_graph(edges)
    # Mark current node as visited and adds to recursion stack
    visited[v] = True
    recStack[v] = True
    # Recur for all neighbours, if any neighbour is visited and in recStack then graph is cyclic
    for neighbour in graph[v]:
        if visited[neighbour] == False:
            if isCyclicUtil(edges, neighbour, visited, recStack) == True:
                return True
        elif recStack[neighbour] == True:
            return True
    # The node needs to be popped from recursion stack before function ends
    recStack[v] = False
    return False


def isCyclic(edges, data):
    """
    Subroutine 2: Test if a given graph is cyclic
    Function is modified from http://www.geeksforgeeks.org/detect-cycle-in-a-graph/
    :param graph: a directed graph in the form of a dictionary
    :return: True if graph is cyclic else false
    """
    numNodes = len(data.columns)
    visited = [False] * numNodes
    recStack = [False] * numNodes
    for node in range(numNodes):
        if visited[node] == False:
            if isCyclicUtil(edges, node, visited, recStack) == True:
                return True
    return False


def obtain_parents(node, edges):
    """
    Get the parent nodes of a specific node
    :param node: Node being evaluated
    :param edges: graph in the form of a list of edges
    :return: a list of the parent nodes
    """
    graph = edges_to_graph(edges)
    parents = []
    for parent, child in graph.items():
        if node in child:
            parents.append(parent)
    return parents


def calculate_cpt(edges, data):
    """
    Subroutine 1: Update the conditional probability table (CPT) for each node given a graph and the data using ML procedure.
    :param edges: graph in the form of a list of edges
    :param data: training data used to build conditional probability table
    :return: two dictionaries, one dictionary holds the cpt of each node in form of a dataframe, the other dictionary
            holds the cpt of each node in the form of a dictIonary (USEFUL FOR DIFFERENT THINGS)
    """
    CPTdict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: -1)))
    CPTdataframes = defaultdict(lambda: -1)
    for node in range(len(data.columns)):
        parents = obtain_parents(node, edges)
        CPTdict[node]['parents'] = parents
        if parents:
            if len(parents) == 1:
                possibilities = [0,1]
            else:
                possibilities = list(itertools.product([0, 1], repeat=len(parents)))
            cpt = pd.DataFrame(possibilities, columns=parents)
            PrB = pd.DataFrame({'PrB': data.groupby(parents).size()}).reset_index()  # P(B)
            PrAB0 = pd.DataFrame({'PrAB0': data.loc[data[node] == 0].groupby(
                parents).size()}).reset_index()  # P(A,B) x=0
            PrAB1 = pd.DataFrame(
                {'PrAB1': data.loc[data[node] == 1].groupby(parents).size()}).reset_index()
            cpt = pd.merge(pd.merge(cpt, PrB, on=parents, how='outer'),
                           pd.merge(PrAB0, PrAB1, on=parents, how='outer'), on=parents, how='outer')
            cpt.fillna(0, inplace=True)
            cpt['Pr0'] = cpt['PrAB0']/cpt['PrB']
            cpt['Pr1'] = cpt['PrAB1']/cpt['PrB']
            cpt.drop(['PrB', 'PrAB0', 'PrAB1'], axis=1, inplace=True)
            index = 0
            for possibility in possibilities:
                CPTdict[node]['Pr0'][possibility] = cpt['Pr0'][index]
                CPTdict[node]['Pr1'][possibility] = cpt['Pr1'][index]
                index += 1
        else:
            cpt = pd.DataFrame(list(itertools.product([0, 1], repeat=len(parents))), columns=parents)
            cpt['Pr0'] = len(data[data[node] == 0]) / len(data)
            cpt['Pr1'] = len(data[data[node] == 1]) / len(data)
            CPTdict[node]['Pr0'] = len(data[data[node] == 0]) / len(data)
            CPTdict[node]['Pr1'] = len(data[data[node] == 1]) / len(data)
        CPTdataframes[node] = cpt
    return CPTdict, CPTdataframes


def propose_new_graph(edges, data):
    """
    Make a random change to a graph, either add an edge, delete an edge, or reverse an edge
    :param edges: a graph in the form of a list of edges
    :param data: training data
    :return: a new graph in the form of a list of edges
    """
    numNodes = len(data.columns)
    proposed_edges = edges[:]
    roulette = random.randint(0,2)
    if roulette == 0:
        # add an edge
        p = random.choice(range(numNodes))
        otherNodes = list(range(0, p)) + list(range(p+1, numNodes))
        c = random.choice(otherNodes)
        proposed_edges.append((p,c))
    if roulette == 1:
        # delete an edge
        edge = random.randint(0,len(edges)-1)
        del proposed_edges[edge]
    if roulette == 2:
        # reverse an edge
        edge = random.randint(0,len(proposed_edges)-1)
        proposed_edges[edge] = (proposed_edges[edge][1], proposed_edges[edge][0])
    proposed_edges.sort()
    return proposed_edges



def isValid(maxParents, edges, data):
    """
    Determine if a proposed graph is valid in that the parents of any nodes don't exceed the maximum parents and checks
    if there are any redundant edges
    :param maxParents: maximum number of parents any node can have
    :param edges: a graph in the form of a list of edges
    :param data: training data
    :return: True if there are no redundant edges and max parents within range, else False
    """
    k = []
    for node in range(len(data.columns)):
        k.append(len(obtain_parents(node, edges)))
    # is the maximum number of parents still 3?
    # is change redundant?
    if max(k) <= maxParents and len(edges) == len(set(edges)):
        return True
    else: return False


def calculate_score(edges, data):
    """
    Calculate the score of a graph
    :param edges: a graph in the form of a list of edges
    :param data: training data
    :return: a negative float (will be negative number in the hundreds, e.g. -350)
    """
    cpt, cptDataframes = calculate_cpt(edges, data)
    score = 0
    for experiment in range(len(data)):
        for node in range(len(data.columns)):
            parents = obtain_parents(node, edges)
            if len(parents) == 1:
                condition = int(data.loc[experiment][parents])
                prob = [cpt[node]['Pr0'][condition] if data.loc[experiment][node] == 0 else cpt[node]['Pr1'][condition]][0]
            elif len(parents) > 1:
                condition = tuple(data.loc[experiment][parents])
                prob = [cpt[node]['Pr0'][condition] if data.loc[experiment][node] == 0 else cpt[node]['Pr1'][condition]][0]
            else:
                prob = [cpt[node]['Pr0'] if data.loc[experiment][node] == 0 else cpt[node]['Pr1']][0]
            score += np.log10(prob)
    return score


def optimize_graph(data, maxParents=3, maxAttempts=50):
    """
    Initialize a graph and make random changes until stopping criteria is met
    :param data: training data
    :param maxParents: maximum number of parents a node can have
    :param maxAttempts: maximum number of valid proposed new graphs that will be evaluated for superiority to current best before stopping
    :return: a dictionary of results, including the starting graph's edges (G0), the starting graph's score (score_G0),
            the best graph's edges (G_best), the best graph's score (score_Gbest), and the best graph's cpt in two forms
            (cpt_Gbest) for use in predict() and (cptGbestDataframes) for printing in a pretty way for homework report,
    """
    print("\nGenerating initial graph...")
    G0 = random_directed_graph(data, maxParents)
    while isCyclic(G0, data):
        G0= random_directed_graph(data, maxParents)
    score_G0 = calculate_score(G0, data)
    print("Initial score: ", score_G0)
    print("Initial graph: ", G0, "\n")
    Gbest = G0[:]
    score_Gbest = score_G0
    while True:
        count = 0
        while True:
            proposed_G = propose_new_graph(Gbest, data)
            while (isCyclic(proposed_G, data)==True) or (isValid(maxParents, proposed_G, data) == False):
                proposed_G = propose_new_graph(Gbest, data)
            proposed_score = calculate_score(proposed_G, data)
            count += 1
            if proposed_score > score_Gbest:
                Gbest = proposed_G[:]
                score_Gbest = proposed_score
                print("Updated score: ",score_Gbest,"\t\tAttempts: ", count)
                print("Updated graph: ", Gbest, "\n")
                count = 0
            if count == maxAttempts:
                break
        if count == maxAttempts:
            print("Final score: ", score_Gbest, "\t\tAttempts: ", count)
            print("Final graph: ", Gbest, "\n")
            cpt_Gbest, cptGbestDataframes = calculate_cpt(Gbest, data)
            break
    return {'G0': G0, 'score_G0': score_G0, 'Gbest': Gbest, 'score_Gbest': score_Gbest, 'cpt_Gbest':cpt_Gbest, 'cptGbestDataframes': cptGbestDataframes}


# - - - - - P A R T   2 - - - - -
# Test and Evaluation:  Test the network by predicting if the expression level of the fifth gene is up (1) or
#                       down (0) given the expression levels of the other genes.


def predict(targetNode, edges, cpt, data):
    """
    Predict the value of a gene given a graph and test data
    :param targetNode: the node you want to predict to determine your graph's prediction accuracy
    :param edges: a graph in the form of a list of edges
    :param cpt: a dictionary of conditional probability tables for the graph being evaluated
    :param data: test data
    :return: a float between 0 and 1 indicating the % prediction accuracy of the gene being evaluated
    """
    parents = obtain_parents(targetNode, edges)
    correct_predictions = 0
    for row in range(len(data)):
        if len(parents) == 1:
            condition = int(data.loc[row][parents])
            prob1 = [cpt[targetNode]['Pr0'][condition] if data.loc[row][targetNode] == 1 else cpt[targetNode]['Pr1'][condition]][0]
            prediction = [1 if prob1 > 0.5 else 0][0]
        if len(parents) > 1:
            condition = tuple(data.loc[row][parents])
            prob1 = [cpt[targetNode]['Pr0'][condition] if data.loc[row][targetNode] == 1 else cpt[targetNode]['Pr1'][condition]][0]
            prediction = [1 if prob1 > 0.5 else 0][0]
        if len(parents) == 0:
            prediction = [1 if cpt[targetNode]['Pr1'] > 0.5 else 0][0]
        if prediction == data[targetNode][row]:
            correct_predictions += 1
    prediction_accuracy = correct_predictions/len(data)
    return prediction_accuracy


def main(trainingDatafile, testDatafile, maxParents, maxAttempts, testGene):
    """
    Main function to generate graph, optimize graph, and analyze prediction accuracy of best graph
    :param trainingDatafile: TSV file of training data in same format as provided for assignment
    :param testDatafile: TSV file of test data in same format as provided for assignment
    :param maxParents: maximum number of parents each node can have
    :param maxAttempts: maximum number of valid proposed new graphs that will be evaluated for superiority to current best before stopping program
    :param testGene: the gene you want to predict to determine your graph's prediction accuracy
    :return:
    """
    train = read_data(trainingDatafile)
    test = read_data(testDatafile)
    results = optimize_graph(train, maxParents=maxParents, maxAttempts=maxAttempts)
    prediction_accuracy = predict(testGene, results['Gbest'], results['cpt_Gbest'], test)
    print("Prediction accuracy for gene %s: " % str(testGene), prediction_accuracy, "\n")
    print("Conditional probability tables of final graph: \n")
    for node in results['cptGbestDataframes']:
        print("Gene: ", node, "\n", results['cptGbestDataframes'][node], "\n")
if __name__ == "__main__":
    main(trainingData, testData, maxParents, maxAttempts, testGene)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - E n d   o f   F i l e - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

