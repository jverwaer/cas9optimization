# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 09:45:05 2017

@author: Jan Verwaeren

@title: module for optimizing the cas9 protein
"""
from collections import Counter
from copy import deepcopy
import numpy as np
from random import shuffle
from itertools import product
import re
import random
###############################################################################
# convenience functions
###############################################################################
def _readHighlyExpressedCDS(fname, check = True):
    
    """ 
    Function to read a file containing several coding domains (one CD per line) and return them as a capitalized list. Line starting with '#' are treated as comments
    
    """
    
    with open(fname, mode = 'r') as f:
        codingDomains_lst = [line.strip().upper() for line in f if len(line) > 5 and line[0] != '#']

    if check:
        for i, seq in enumerate(codingDomains_lst):
            if len(seq) % 3 != 0 or  (seq.count("A") + seq.count("C") + seq.count("T") + seq.count("G")) != len(seq):
                raise ValueError("Sequence {:d} is not valid!".format(i+1))
        
    return codingDomains_lst

def _createCodonTable():
    """
    Function to create a dictionaur in which the keys are codons and the values are amino acids
    """
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table

def _createAATable():
    """
    Function to create a dictionary in which keys ara amino acids and values are a list of codons (the stop codon is translated as '*')
    """
    AA_TO_CODON = dict()
    for codon, AA in CODON_TO_AA.items():
        try:
            AA_TO_CODON[AA].append(codon)
        except KeyError:
            AA_TO_CODON[AA] = [codon]
        
    return AA_TO_CODON


def DNASeq2CodonList(dnaSeq):
    """
    Function to tranlate a DNA string to a list of codons
    """
    return [dnaSeq[i:i+3] for i in range(0,len(dnaSeq),3)]

def translateDNA(dnaSeq):
    """
    Function translates a DNA sequence (string or list of codons) into a AA sequence (returns a string). The last character is a '*' which codes for the stop codon
    """
    if type(dnaSeq) is str:
        AAseq = [CODON_TO_AA[dnaSeq[i:i+3]] for i in range(0,len(dnaSeq),3)]
        
    elif type(dnaSeq) in [list, tuple] and len(dnaSeq[0]) == 3:
        AAseq = [CODON_TO_AA[codon] for codon in dnaSeq]
        
    else:
        raise ValueError("Unsupported input for 'translateDNA', should eather be a DNA string or a list of codons")
    
    return("".join(AAseq))


def reverseAAseq(AAseq):
    """
    Function re-translates an AA sequence into a list of codon-lists
    """
    return [AA_TO_CODON[AA] for AA in AAseq]
    
        
def approximatePMF(targetPMF, sampleSize):
    """
    Function computes a frequency distribution with specified sample size (np.sum(result) == sampleSize) zich that it respects the targetPMF as closely as possible
    """
    targetPMF = targetPMF/np.sum(targetPMF) * sampleSize
    samplePMF = np.array(np.floor(targetPMF), dtype=int)
    shortage = int(sampleSize - np.sum(samplePMF))
    for i in range(shortage):
        samplePMF[np.argmax(targetPMF - samplePMF)] += 1
    return samplePMF
    
def _createAApiarTable():
    aaPairToCodon = dict()
    for AA1, codon1 in AA_TO_CODON.items():
        for AA2, codon2 in AA_TO_CODON.items():
            aaPairToCodon[AA1+AA2] = [element[0]+element[1] for element in list(product(codon1, codon2))]
    return aaPairToCodon

def readCodonFrequenciesFromFile(fname = "codon_usage.txt"):
    """
    Specific function for reading codon frequency table from a file like http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=36050
    Not to be used frequently
    
    """
    with open(fname, mode = 'r') as f:
        table = f.read()
    
    table= table.replace("U", "T")
    table= table.split("\n")
    
    codonFreq = []
    for i in range(len(table)): 
        if len(table[i]) > 5:
            codonFreq.append((table[i][0:3], float(table[i][3:8])))
            codonFreq.append((table[i][18:21], float(table[i][21:26])))
            codonFreq.append((table[i][36:39], float(table[i][39:44])))
            codonFreq.append((table[i][54:57], float(table[i][57:62])))
    codonFreq = dict(codonFreq)
    return CodonFrequency({ AA: np.array([ codonFreq[codon] for codon in AA_TO_CODON[AA]])  for AA in AA_TO_CODON})

def findOccurences(mystr, templates_lst):
    """
    Function to find all occurences of the strings in templates_lst in mystr. This function returns a sorted list of tuples (index, string).
    """
    
    indices = []
    for template in templates_lst:
        regExp = '(?='+template+')'
        indices_templ = [m.start() for m in re.finditer(regExp, mystr)]
        indices.extend([(index, template) for index in indices_templ])
    
    indices.sort()
    return indices
    
def frequencyBasedSort(values_lst, freq_lst):
    freq_values= sorted([(freq_lst[i], value) for i, value in enumerate(values_lst)], reverse = True)
    return [freq_val[1] for freq_val in freq_values]

###############################################################################
# convenience variables (constants)
###############################################################################
CODON_TO_AA = _createCodonTable()
AA_TO_CODON = _createAATable()
AAPAIR_TO_CODON = _createAApiarTable()


###############################################################################
# class defenitions
###############################################################################
class CodingDomain(object):
    """
    Class represents a coding domain --internally stored as a list of codons. Takdes a dna-string as input of the initialiser, or a list of codons.
    
    """
    def __init__(self, dnaSeq, computeCodonFreq = True):
        
        """
        PARAMETERS
        ---------
        dnaSeq : str or list
            
            a DNA sequence or a list of codons
        
        """
        
        if type(dnaSeq) is list and len(dnaSeq[0]) == 3:
            
            self.__codonlist = dnaSeq
        else:
            self.__codonlist = DNASeq2CodonList(dnaSeq)
            
        if computeCodonFreq:
            self.computeCodonFrequency()
            self.computeCodonPairFrequency()
        
    def computeCodonFrequency(self):
        self.__codonFrequency_Counter = Counter(self.__codonlist)
        self.__codonFrequency = CodonFrequency([(AA , (np.array([self.__codonFrequency_Counter[codon] for codon in codonList]))) for AA, codonList in AA_TO_CODON.items()])
    
    def get_codonFrequency(self, as_counter = False):
        
        if as_counter:
            return self.__codonFrequency_Counter
        
        return self.__codonFrequency
    
    def get_codonList(self):
        return self.__codonlist
    
    def as_AAstring(self):
        return translateDNA(self.__codonlist)
    
    def as_AAsequence(self):
        return AAsequence(translateDNA(self.__codonlist))
    
    def as_DNAstring(self):
        return "".join(self.__codonlist)
    
    def __str__(self):
        return self.__codonlist.__str__()
    
    def computeCodonPairFrequency(self):
        self.__codonPairFrequency_Counter = Counter(["".join(self.__codonlist[i:i+2]) for i in range(len(self.__codonlist)-1)])
        self.__codonPairFrequency = CodonPairFrequency([(AApair , (np.array([self.__codonPairFrequency_Counter[codon] for codon in codonList]))) for AApair, codonList in AAPAIR_TO_CODON.items()])
        
    def get_codonPairFrequency(self, as_counter = False):
        if as_counter:
            return self.__codonPairFrequency_Counter
        
        return self.__codonPairFrequency
    
    def changeCodon(self, index, newCodon, sameAA = True):
        """
        Method to change the codon at position index (an int) into newCodon (a string). If sameAA is true, a check is performed to see if the codon encodes the same AA. This method also UPDATES the codonfrequency atrributes.
        
        RETURNS
        ------
        List of (at most) three tuples, containing the AA, and its new codonFrequency, the remaining two contain the same info for the pairs codonPairFrequencies
        
        """
        oldCodon = self.__codonlist[index]
        AA = CODON_TO_AA[oldCodon]
        candidateCodons = AA_TO_CODON[AA]
        
        returnList = []
        
        if sameAA:
            try:
                # change codonFrequency
                #store_old_AA = np.copy(self.__codonFrequency[AA])
                self.__codonFrequency[AA][candidateCodons.index(oldCodon)] -= 1
                self.__codonFrequency[AA][candidateCodons.index(newCodon)] += 1
                returnList.append((AA, self.__codonFrequency[AA]))
                # change codonPairFrequency
                if index > 0:
                    min1Codon = self.__codonlist[index-1]
                    min1AA = CODON_TO_AA[min1Codon]
                    candidatePairs_min1 = AAPAIR_TO_CODON[min1AA+AA]
                    #store_old_min1AA_AA = np.copy(self.__codonPairFrequency[min1AA+AA])
                    self.__codonPairFrequency[min1AA+AA][candidatePairs_min1.index(min1Codon+oldCodon)] -= 1
                    self.__codonPairFrequency[min1AA+AA][candidatePairs_min1.index(min1Codon+newCodon)] += 1
                    returnList.append((min1AA+AA, self.__codonPairFrequency[min1AA+AA]))
                else:
                    returnList.append((0,))
                    
                if index < (len(self.__codonlist) - 1):
                    plus1Codon = self.__codonlist[index+1]
                    plus1AA = CODON_TO_AA[plus1Codon]
                    candidatePairs_plus1 = AAPAIR_TO_CODON[AA+plus1AA]
                    #store_old_AA_plusAA= np.copy(self.__codonPairFrequency[AA+plus1AA])
                    self.__codonPairFrequency[AA+plus1AA][candidatePairs_plus1.index(oldCodon+plus1Codon)] -= 1
                    self.__codonPairFrequency[AA+plus1AA][candidatePairs_plus1.index(newCodon+plus1Codon)] += 1
                    returnList.append((AA+plus1AA , self.__codonPairFrequency[AA+plus1AA]))
                else:
                    returnList.append((0,))
                # make update
                self.__codonlist[index] = newCodon
                return returnList # returns a list containing the ols and new AA frequencies
            except ValueError:
                raise ValueError("Probably, the new codon codes for a different AA than the original, which is not allowad as sameAA == True.")
        else:
            raise NotImplementedError("Only changes that do not affect the AA sequence have been implemented thus far")
        
    def hammingDistance(self, cd):
        """
        Method to compute hamming distance between two coding domains
        """
        hamming_distance = 0
        for i in range(len(self.get_codonList())):
            if self.get_codonList()[i] != cd.get_codonList()[i]:
                hamming_distance += 1
        return hamming_distance
                
    def replaceCodons(self,replacelist):
        """
        Methods that uses the ('old_codon', 'new_codon') pairs in replacelist such that all occurences of 'old_codon' are replaced with 'new_codon'. 
        
        PARAMETERS
        ----------
        replacelist : list
            List of the format [('TTG', 'CTC'), ('ACA', 'ACC'), ...] consisting of old_codon (index 0) en new_codon (index 1) tuples.
        
        """
        
        replaceDict = dict(replacelist)
        
        for i, codon in enumerate(self.__codonlist):
            if codon in replaceDict:
                self.changeCodon(i, replaceDict[codon])
        
# ----------------------------------------------------------------------------------------------------------------------------------------

class CodonFrequency(dict):
    """
    Class that implements a codon frequency. This is simply a dict (with some additional methods). Keys are AA's and the values a list of codon frequencies. Use the getters and setters for interacting with this class. 
    """
    def __init__(self, *args, **kwargs):
        
        if len(args) > 0 and type(args[0]) is Counter:
            counter = args[0]
            self.__init__([(AA , (np.array([counter[codon] for codon in codonList]))) for AA, codonList in AA_TO_CODON.items()])
        else:
            super(CodonFrequency, self).__init__(*args, **kwargs)
        
    def as_codonFrequency(counter):
        return CodonFrequency([(AA , (np.array([counter[codon] for codon in codonList]))) for AA, codonList in AA_TO_CODON.items()])
    
    def __sub__(self, freq):
        return CodonFrequency([(AA , np.abs(self[AA] - freq[AA])) for AA in AA_TO_CODON])
    
    def __abs__(self):
        freq_sum = 0
        for AA, freq in self.items():
            freq_sum += np.sum(freq)
        return freq_sum
    
    def normalize(self):
        """
        Normalizes this object
        """
        for AA in self:
            if np.sum(self[AA]) > 0.0001:
                self[AA] = self[AA] / np.sum(self[AA])
            else:
                self[AA] = np.ones(len(self[AA])) / len(self[AA])
    
    def normalizedVersion(self):
        """
        Returns normalized version, original object remains unchanged
        """
        c = deepcopy(self)
        c.normalize()
        return c
    
    
    def lowFrequentCodonsAndReplacements(self, threshold, normalized = False):
        
        lowFreqCodons_lst = []
        
        if normalized:
            tmp_codonFreq = self.normalizedVersion()
        else:
            tmp_codonFreq = self
        
        for AA, codonlist in AA_TO_CODON.items():
            
            for i, aaFreq in enumerate(tmp_codonFreq[AA]):
                if aaFreq < threshold:
                    lowFreqCodons_lst.append((codonlist[i], codonlist[np.argmax(tmp_codonFreq[AA])]))
                    
        return lowFreqCodons_lst

# ----------------------------------------------------------------------------------------------------------------------------------------

class CodonPairFrequency(dict):
    
    def __init__(self, *args, **kwargs):
        
        if len(args) > 0 and type(args[0]) is Counter:
            counter = args[0]
            self.__init__([(AApair , (np.array([counter[codon] for codon in codonList]))) for AApair, codonList in AAPAIR_TO_CODON.items()])
        else:
            super(CodonPairFrequency, self).__init__(*args, **kwargs)
            
    def __sub__(self, freq):
        return CodonPairFrequency([(AApair , np.abs(self[AApair] - freq[AApair])) for AApair in AAPAIR_TO_CODON])
    
    def __abs__(self):
        freq_sum = 0
        for AApair, freq in self.items():
            freq_sum += np.sum(freq)
        return freq_sum
    
    def normalize(self):
        for AApair in self:
            if np.sum(self[AApair]) > 0.0001:
                self[AApair] = self[AApair] / np.sum(self[AApair])
            else:
                self[AApair] = np.ones(len(self[AApair])) / len(self[AApair]) 

# ----------------------------------------------------------------------------------------------------------------------------------------

class CDcreationInfo(dict):
    
    """
    Class that contains info on how to derive a coding domain from a given AA-sequence
    """
    
    def __init__(self, method = "random", codonFrequency = None, sumToOneCheck = True):
        """
        PARAMETERS
        ---------
        method : str
            should be one of 'random', 'codonfrequency_soft' or 'codonfrequency_hard'
            
        codonFrequency : CodonFrequency
            should be a codonfrequency or None (ignored if method = 'random')
            
        sumToOneCheck : bool
            True if the codonFrequency should be transformed into probabilities, False if not (use false only in case of efficiency issues)
        
            
        """
        
        dict.__init__(self, [('method', method), ('codonFrequency', codonFrequency), ('sumToOneCheck', sumToOneCheck)])


DEFAULT_INFO = CDcreationInfo()

# ----------------------------------------------------------------------------------------------------------------------------------------

class AAsequence(object):
    """
        Class that implements Amino Acid sequences. Als string methods can be used.
        
        ATTRIBUTES
        ----------
        
        sequence : str
        
            The string containing the AA sequence
        
        codingDomain : CodingDomain
        
            A coding domain for the AA sequence
    """
    
    def __init__(self, sequence):
        """
        PRAMETERS
        ---------
        
        sequence : any sequence of AA-characters that can be cast into a string (typically just a string of AA's)
        
        """
        self.sequence = str(sequence)
        
    
    def set_codingDomain(self, codingDomain):
        self.__codingDomain = codingDomain
        
    def get_codingDomain(self):
        try:
            return self.__codingDomain
        except AttributeError:
            raise AttributeError('The coding domain has not be set yet, use create_codingDomain or set_codingDomain first! ')
    
    def create_codingDomain(self, createInfo =  DEFAULT_INFO , addAsAttribute = True):
        """
            Method to generate a coding domain for this AA sequence.
            
            PARAMETERS
            ----------
            createInfo : dict
                dict containing the info to generate a coding domain for this AA sequence. 
            
        """
        if createInfo['method'] == 'random':
            self.__codingDomain = self.__codingDomainFromCodonFrequency(codonFrequency = "uniform")
        elif createInfo['method'] == 'codonfrequency_soft':
            self.__codingDomain = self.__codingDomainFromCodonFrequency(codonFrequency = createInfo['codonFrequency'], independentsampling=True, sumToOneCheck = createInfo['sumToOneCheck'])
        elif createInfo['method'] == 'codonfrequency_hard':
            self.__codingDomain = self.__codingDomainFromCodonFrequency(codonFrequency = createInfo['codonFrequency'], independentsampling=False, sumToOneCheck = createInfo['sumToOneCheck'])
            
            
    def __codingDomainFromCodonFrequency(self, codonFrequency = "uniform", independentsampling = True, sumToOneCheck = True):
        
        if codonFrequency == "uniform": # sample uniformly from available codons
            codingDomain = CodingDomain([np.random.choice(AA_TO_CODON[AA]) for AA in self.sequence])
            
        else: # use distribution provided in codonFrequency
        
            if sumToOneCheck:
                
                codonFrequency = {AA : codonFrequency[AA]/np.sum(codonFrequency[AA]) for AA in codonFrequency if np.sum(codonFrequency[AA]) > 0}  # transform frequencies into probabilities
        
            if independentsampling: # independent sample according to PMF
                
                codingDomain = CodingDomain( [np.random.choice(AA_TO_CODON[AA], p = codonFrequency[AA]) for AA in self.sequence])
                
            else: # resulting sample should match target pmf as closely as possible
                AAfrequency = dict(Counter(self.sequence))
                AAfrequency = {AA : list(np.repeat(AA_TO_CODON[AA], approximatePMF(codonFrequency[AA], AAfrequency[AA]))) for AA in AAfrequency}
                for AA, codonlist in AAfrequency.items():
                    shuffle(codonlist)
                codingDomain = CodingDomain ([AAfrequency[AA].pop() for AA in self.sequence])
        
        return codingDomain
    
    def get_codonFrequency(self):
        return self.get_codingDomain().get_codonFrequency()






# ----------------------------------------------------------------------------------------------------------------------------------------       

def unittest():
    # create simple sequence
    #Eseq = AAsequence("ACDEFCDAAFFE")
    seq = AAsequence('L'*1000 + 'N' * 1000)
    seq.create_codingDomain()
    assert seq.get_codingDomain().as_AAstring() == seq.sequence, "Uniform sampling failed"
    
    import os
    os.chdir("E:\\internAdvies\\jiang_tan")
    extractor = CodonFrequencyExtractor('highly_expressed_cds.txt')
    
    creationInfo = CDcreationInfo(method = 'codonfrequency_hard', codonFrequency=extractor.codonFrequency)
    
    seq.create_codingDomain(creationInfo)
    assert seq.get_codingDomain().as_AAstring() == seq.sequence, "Frequency-based sampling failed"
    
    return seq.get_codonFrequency()
    
    
# ----------------------------------------------------------------------------------------------------------------------------------------       
            

class CodonFrequencyExtractor(object):
    """
        Class that implements a frequency extractor of codons (frequency distribution for of each codon per AA is extracted from a list of highly expressed coding domains)
        
        ATTRIBUTES
        ----------
        codingDomains : list of strings
            List of strings contining the highly expressed coding domains 
        
        codonFrequency : dict (codon - absolute frequency (int))
            Dict of codon frequencies
            
        codonPairFrequency : dict (codonpairs - absolute frequency (int))
            Dict of codon pair frequencies
        
    """
    def __init__(self, codingDomains):
        """
            
            PARAMETERS
            ----------
            codingDomains : str (filename) or list of strings
                List of strings containing the coding domains or a filename referring to a file that contains coding domains
        """
        if type(codingDomains) is str:
            codingDomains = _readHighlyExpressedCDS(codingDomains)
            
        self.codingDomains  = [CodingDomain(dnaSeq) for dnaSeq in codingDomains]
        self.computeFrequencyCondons()
        self.computeFrequencyCodonPairs()
    
    def readHighlyExpressedCDS(self, fname):
        codingDomains = _readHighlyExpressedCDS(fname)
        self.codingDomains  = [CodingDomain(dnaSeq) for dnaSeq in codingDomains]
        self.computeFrequencyCondons()

    def computeFrequencyCondons(self):
        som = Counter()
        for cdom in self.codingDomains:
            som = som + cdom.get_codonFrequency(True)
        # compute absolute frequencies per AA
        self.codonFrequency = CodonFrequency.as_codonFrequency(som)
        
    def computeFrequencyCodonPairs(self):
        som = Counter()
        for cdom in self.codingDomains:
            som = som + cdom.get_codonPairFrequency(True)
        # compute absolute frequencies per AA
        self.codonPairFrequency = CodonPairFrequency(som)
       
# ----------------------------------------------------------------------------------------------------------------------------------------       
       

class CodonFrequencyAnalyser(CodonFrequencyExtractor):
    
    """
    Sublcass of CodonFrequencyExtractor, with the purpose of building an additional distance matrix that captures the difference between the codon frequencies for the different coding domains.
    
    ONLY TO BE USED FOR COMPAIRING (PAIR) FREQUENCIES --> as frequencies are normalised
    
    ATTRIBUTES
    ----------
    
    codingDomains : list of strings
            List of strings contining the highly expressed coding domains 
        
    codonFrequency : dict (codon - absolute frequency (int))
            Dict of codon frequencies
    
    frequencyDistanceMatrix : np.array
            Matrix of distances between codon frequencies in different coding domains 
    
    frequencyPairsDistanceMatrix : np.array
            Matrix of distances between codonpair frequencies in different coding domains 
    
    """
    
    def __init__(self, codingDomains):
        
        CodonFrequencyExtractor.__init__(self, codingDomains)
        
        self.frequencyDistanceMatrix = np.zeros((len(self.codingDomains), len(self.codingDomains)))
        self.frequencyPairsDistanceMatrix = np.zeros((len(self.codingDomains), len(self.codingDomains)))
        
        for codingDomain in self.codingDomains:
            codingDomain.get_codonFrequency().normalize()
            codingDomain.get_codonPairFrequency().normalize()
        
        for i, codingDomain1 in enumerate(self.codingDomains):
            
            for j, codingDomain2 in enumerate(self.codingDomains):
                
                self.frequencyDistanceMatrix[i,j] = abs(codingDomain1.get_codonFrequency() - codingDomain2.get_codonFrequency())
                self.frequencyPairsDistanceMatrix[i,j] = abs(codingDomain1.get_codonPairFrequency() - codingDomain2.get_codonPairFrequency())
                
                
    
    
# ----------------------------------------------------------------------------------------------------------------------------------------       
    
class CodonOptimizer(object):
    """
    Class for finding an optiimal coding domain for a given AA sequence.
    """
    
    def __init__(self, aaSequence, codonFrequency_target, codonPairFrequency_target):
        """
        PARAMETERS
        ----------
        
        aaSequence : str or AAsequence
            the AA sequence of the protein for which an optimal coding domain is needed
            
        codonFrequency_target : CondonFrequency
            Codon frequency observed for the organism
            
        codonPairFrequency_target : CodonPairFrequency
            Frequency of the codon pairs for the target
        
        """
        if type(aaSequence) is AAsequence: # ensures that a copy is made
            self.aaSequence = AAsequence(aaSequence.sequence)
        else:
            self.aaSequence = AAsequence(aaSequence)
        
        self.codonFrequency_target, self.codonPairFrequency_target = self.rescaleAccordingToAAsequence(codonFrequency_target, codonPairFrequency_target)
        self.nAA = len(self.aaSequence.sequence)
    
    def rescaleAccordingToAAsequence(self, codonFrequency_target, codonPairFrequency_target):
        """
        Method to rescale the codonfrequencies (and pairs) according to the current AA sequence
        """
        
        rescaledCodonFrequency = dict()
        for AA in AA_TO_CODON:
            count_target = np.sum(codonFrequency_target[AA]) # count nr of occurence of AA in target 
            count_current = np.sum(self.aaSequence.sequence.count(AA)) # count nr of occurence of AA in current 
            if count_target != 0:
                rescaledCodonFrequency[AA] = codonFrequency_target[AA] / count_target * count_current
            else:
                rescaledCodonFrequency[AA] = np.ones(len(codonFrequency_target[AA])) / len(codonFrequency_target[AA]) * count_current
                
        rescaledCodonPairFrequency = dict()
        for AApair in AAPAIR_TO_CODON:
            count_target = np.sum(codonPairFrequency_target[AApair])
            count_current = len(re.findall('(?={:s})'.format(AApair.replace('*', '\*')), self.aaSequence.sequence))
            if count_target != 0:
                rescaledCodonPairFrequency[AApair] = codonPairFrequency_target[AApair] / count_target * count_current
            else:
                rescaledCodonPairFrequency[AApair] = np.ones(len(codonPairFrequency_target[AApair])) / len(codonPairFrequency_target[AApair]) * count_current
        
        return rescaledCodonFrequency, rescaledCodonPairFrequency
        
    def optimizeCodingDomainUnivariate(self, createInfo = None):
        
        if createInfo is None:
            createInfo = CDcreationInfo(method = "codonfrequency_hard", codonFrequency = self.codonFrequency_target)
        
        self.aaSequence.create_codingDomain(createInfo =  createInfo, addAsAttribute = True)
        
    def AAatIndex(self, index):
        return self.aaSequence.sequence[index]
    
    def optimizeCodingDomainBivariate(self, initialGuess = 'univariate', nIter = 100, algorithm = "HC", annealingTemperatures = list(np.linspace(1, 0, 8))*10):
        """
        Function to find optimal coding domain, using second order info.
        
        PARAMETERS
        ----------
        
        initialGuess : 'univariate', 'random' or a CodingDomain
        
        nIter : int or None
            Number of iterations, if None, the procedure iterates untill convergence
            
        algorithm : str
            The string 'HC' for hill climing or 'SA' for simulate annealing (in case of SA, also indicate the annealing temperatures)
        """
        
        if initialGuess == "univariate":
            self.optimizeCodingDomainUnivariate()
            
        elif initialGuess == 'random':
            createInfo = CDcreationInfo(method = "random")
            self.optimizeCodingDomainUnivariate(createInfo = createInfo)
            
        elif type(initialGuess) is CodingDomain:
            self.aaSequence.set_codingDomain(initialGuess)
        
        else:
            raise ValueError("Invalid choice for the initial guess")
        
        self.initialGuess = deepcopy(self.aaSequence.get_codingDomain().get_codonPairFrequency())
        
        if algorithm == "HC":
            message = self.hillClimbing(nIter)
        elif algorithm == "SA":
            message = self.simulatedAnnealing(annealingTemperatures = annealingTemperatures)
        else:
            raise ValueError("Algorithm not specified correctly")
        
        return message
    
    def hillClimbing(self, nIter):
        """
        Hill climing optimizer
        """
        iteration = 0
        converged = False
        while iteration < nIter and (not converged):
            converged = True
            iteration += 1
            for index in range(self.nAA): 
                
                change_applied = self.optimizeStep(index)
    
                if change_applied:
                    converged = False
                    
        message = {"Converged" : converged, "Iterations" : iteration}
        
        return message
    
    
    def simulatedAnnealing(self, annealingTemperatures = list(np.linspace(1.2, 0, 10))*10):
        
        iteration = 0
        for temp in annealingTemperatures:
            iteration += 1
            tolerances = np.random.exponential(scale = temp, size = self.nAA)
            for index in range(self.nAA): 
                
                self.optimizeStep(index, tolerances[index], randomShuffle = True)
                    
        message_hc = self.hillClimbing(100)
        message = {"Converged" : message_hc["Converged"], "Iterations" : iteration}
        return message
    
    def optimizeStep(self, index, tolerance = 0.0, print_info = False, randomShuffle = False):
        
        """
        Method to attempt to optimize a single position (parameter index) in the coding domain. The objective is to find a replacement such that the dissimilarity with codonPairFrequency_target decreases. 
        Which means that: objective_new < objective_old + tolerance (if tolerance == 0, default) this means traditional hill climbing. Once a replacement has is found that respects the former condition, it is applied to the coding domain.
        
        RETURNS
        -------
        
        True in case a change has been applied, False in case no change was applied
        
        """
        
        #AA_freq_target = self.codonFrequency_target[self.AAatIndex(index)]
        #AA_freq_cd = np.copy(self.get_codingDomain().get_codonFrequency()[self.AAatIndex(index)])

        objective = 0.0
        
        if index > 0:
            minAA_AA_pair_freq_target = self.codonPairFrequency_target[self.AAatIndex(index-1)+self.AAatIndex(index)]
            minAA_AA_pair_freq_cd = self.get_codingDomain().get_codonPairFrequency()[self.AAatIndex(index-1)+self.AAatIndex(index)]
            objective += np.sum(np.abs(minAA_AA_pair_freq_target - minAA_AA_pair_freq_cd))
            
        if index < (len(self.aaSequence.sequence) - 1):
            AA_plusAA_pair_freq_target = self.codonPairFrequency_target[self.AAatIndex(index)+self.AAatIndex(index+1)]
            AA_plusAA_pair_freq_cd = self.get_codingDomain().get_codonPairFrequency()[self.AAatIndex(index)+self.AAatIndex(index+1)]
            objective += np.sum(np.abs(AA_plusAA_pair_freq_target - AA_plusAA_pair_freq_cd))
        
        if print_info: # for debugging purposes
            print("START objective:",objective)
            print("Target-freq min:",minAA_AA_pair_freq_target)
            print("CD-freq min:", minAA_AA_pair_freq_cd)
            print("Target-freq plus:",AA_plusAA_pair_freq_target)
            print("CD-freq plus:", AA_plusAA_pair_freq_cd)
        
        if randomShuffle:
            alternative_set = AA_TO_CODON[self.AAatIndex(index)].copy()
            shuffle(alternative_set)
        else:
            alternative_set = AA_TO_CODON[self.AAatIndex(index)]
        
        for alternativeCodon in alternative_set:
            # initial codon
            initial_codon = self.get_codingDomain().get_codonList()[index]
            if initial_codon == alternativeCodon:
                continue
            # apply change and get info of change on frequencies
            change_info = self.get_codingDomain().changeCodon(index, alternativeCodon)
            new_objective = 0
            if index > 0:
                new_objective+=np.sum(np.abs(minAA_AA_pair_freq_target-change_info[1][1]))
            if index < (len(self.aaSequence.sequence) - 1):
                new_objective+=np.sum(np.abs(AA_plusAA_pair_freq_target-change_info[2][1]))
            
            if print_info: # for debugging purposes
                print("-"*20)
                print("NEW objective:",new_objective)
                print("Target-freq min:",minAA_AA_pair_freq_target)
                print("CD-freq min:", change_info[1][1])
                print("Target-freq plus:",AA_plusAA_pair_freq_target)
                print("CD-freq plus:", change_info[2][1])
            
            #obj_diff.append(new_objective-objective)
            if new_objective < (objective - 0.0001) + tolerance: # check if similarity between target and new coding domain is better
                # print("Change applied")
                return True
            else:
                self.get_codingDomain().changeCodon(index, initial_codon)
        return False
        
    def get_codingDomain(self):
        return self.aaSequence.get_codingDomain()
        

    def replaceUnwantedSubstrings(self, unwanted_lst):
        """
        Method that replaces codons such that the dna strings in unwanted_lst are eliminated from the codingDomain
        """
        
        # find ocurrence of substring
        cd_str = self.get_codingDomain().as_DNAstring()
        
        # indices
        indices_unwanted_lst = findOccurences(cd_str, unwanted_lst)
        while len(indices_unwanted_lst) > 0:
            firstIndex, firstTemplate = indices_unwanted_lst[0]
            indices_unwanted_lst = self.__replaceUnwanted_single(firstIndex, firstTemplate, unwanted_lst)
            
                
    def __replaceUnwanted_single(self, firstIndex, template, unwanted_lst):
        
        for codonlistIndex in range(firstIndex//3, ((firstIndex + len(template)-1)//3) + 1):
            original_codon = self.get_codingDomain().get_codonList()[codonlistIndex]
            
            # find candidate replacements (index in __codingDomain and alternatives)
            alternatives = AA_TO_CODON[CODON_TO_AA[self.get_codingDomain().get_codonList()[codonlistIndex]]]
            freq_array = self.codonFrequency_target[CODON_TO_AA[self.get_codingDomain().get_codonList()[codonlistIndex]]]
            
            alternatives = frequencyBasedSort(alternatives, freq_array)
            
            for codon in alternatives:
                # do replace
                self.get_codingDomain().changeCodon(codonlistIndex, codon)
                cd_str = self.get_codingDomain().as_DNAstring()
                indices_unwanted_lst = findOccurences(cd_str, unwanted_lst)
                # check if the replacement solves the problem
                if len(indices_unwanted_lst) == 0 or indices_unwanted_lst[0][0] > firstIndex:
                    print("Position:", codonlistIndex*3, original_codon, "replaced with", codon)
                    return indices_unwanted_lst
                else:
                    self.get_codingDomain().changeCodon(codonlistIndex, original_codon)
        
        raise ReplacementImpossibleError("No valid replacement found")
    
    
#------------------------------------------------------------------------------
class ReplacementImpossibleError(ValueError):
    pass


        