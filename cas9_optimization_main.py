
import cas9_optimization as cas

###############################################################################
# STEP 1: obtain codon frequencies of higly expressed coding domains
#         and optimize (in this example for CAS9 but can be any coding domain)
###############################################################################

# construct a condonfrequency-extractor: the constuctor takes a tekst file
# containing a set of higly expressed coding domains as input (one coding 
# domain per line)
extractor = cas.CodonFrequencyExtractor('data/highly_expressed_cds.txt')
# dictionary ofAAs (keys) and related codons (values)
print(cas.AA_TO_CODON)
# frequencies
print(extractor.codonFrequency)

# read the cas9 coding domain (DNA-form)
with open("data/cas9_coding_domain_T_reesei.txt", mode = "r") as f:
    DNAseq = f.readlines()[1].strip()

# transform into AA sequence
AAseq = cas.translateDNA(DNAseq)

# optimize the codon sequencies (repeat 10 times due to stochasticity in
# the optimization procedure - which is Simulated Annealing) 
results = []
for i in range(10):
    # create an optimizer object
    optimizer = cas.CodonOptimizer(AAseq, extractor.codonFrequency,
                                   extractor.codonPairFrequency)
    # perform codon optimization where the frequency distribution of codon pairs
    # is choosen to match the frequencey observed in the highly expressed domains
    message = optimizer.optimizeCodingDomainBivariate(algorithm="SA")
    
    # compute difference between initialGuess and target
    obj_initial = abs(optimizer.initialGuess - optimizer.codonPairFrequency_target)
    obj_final = abs(optimizer.get_codingDomain().get_codonPairFrequency() -
                    optimizer.codonPairFrequency_target)
    
    results.append((abs(obj_initial),   # initial value of objective function
                    abs(obj_final),     # final value of objective function (low is good)
                    abs(obj_initial) - abs(obj_final), #  improvement due to optimization
                    AAseq == optimizer.get_codingDomain().as_AAstring(), # check if the AA seq has not been changed by accident
                    message["Converged"],
                    message['Iterations'],
                    optimizer.get_codingDomain())) # the coding domain after optimizing
    print("Run:", i+1)

###############################################################################
# STEP 2: do some checks
###############################################################################

# A ) Compute how well the 'cas9_coding_domain_T_reesei' sequence already
# respects the codon frequencies observed in the higly expressed CDs
# ---------------------------------------------------------------------------
    
# create a coding domain object
CD_orig = cas.CodingDomain(DNAseq)
# create an optimizer object
optimizer = cas.CodonOptimizer(AAseq, extractor.codonFrequency ,extractor.codonPairFrequency)

# extract codon frequency and codon pair frequency
codonFrequency_orig = CD_orig.get_codonFrequency()
codonPairFrequency_orig = CD_orig.get_codonPairFrequency()

# (1) compute the difference between the observed codon frequencies in 'cas9_coding_domain_T_reesei'
# and the observed frequencies in the highly expressed coding domains
# (2) identical to (1) but for frequencies of pairs
orig_vs_target_aaFreq = abs(cas.CodonFrequency(optimizer.codonFrequency_target) - codonFrequency_orig)
orig_vs_target_pairsFreq = abs(cas.CodonPairFrequency(optimizer.codonPairFrequency_target) - codonPairFrequency_orig)



# B ) Identical to (A) but with optimized frequences
# ---------------------------------------------------------------------------

# optimized is the 'best' out of the 10 runs in STEP 10
optimized = results[[x[1] for x in results].index(min([x[1] for x in results]))][6]
# compute differences
optimized_vs_target_aaFreq = abs(cas.CodonFrequency(optimizer.codonFrequency_target) - optimized.get_codonFrequency())
optimized_vs_target_pairsFreq = abs(cas.CodonPairFrequency(optimizer.codonPairFrequency_target) - optimized.get_codonPairFrequency())

# print string version of optimized
optimized_str = "".join(optimized.get_codonList())
print(optimized_str)

# PRINT the results
print("SUMMARY RESULTS")
print("Distance between absolute codon pair frequencies observed", 
      "in 'cas9_coding_domain_T_reesei' and frequencies =", orig_vs_target_pairsFreq)



# =============================================================================
# # ALTERNATIVE APPROACH: optimize codon frequencies (AS OPPOSED TO THE CODON
# # PAIR FREQUENCIES) such that the frequencies math those observed in the
# # higly expressed coding domains
# =============================================================================

# optimizer.optimizeCodingDomainUnivariate()
# optimized_UNIVAR = optimizer.get_codingDomain()

# optimized_UNIVAR_vs_target_aaFreq = abs(cas.CodonFrequency(optimizer.codonFrequency_target) - optimized_UNIVAR.get_codonFrequency())
# optimized_UNIVAR_vs_target_pairsFreq = abs(cas.CodonPairFrequency(optimizer.codonPairFrequency_target) - optimized_UNIVAR.get_codonPairFrequency())

# # compute hamming distances
# hamming_distance_UNIVAR_orig = optimized_UNIVAR.hammingDistance(CD_orig)
# hamming_distance_orig = optimized.hammingDistance(CD_orig)
# hamming_distance_alternatives = optimized.hammingDistance(optimized_UNIVAR)

# optimizer2 = CodonOptimizer(AAseq, extractor.codonFrequency ,extractor.codonPairFrequency)
# optimizer2.optimizeCodingDomainUnivariate()
# optimized_UNIVAR2 = optimizer2.get_codingDomain()
# optimized_UNIVAR.hammingDistance(optimized_UNIVAR2)


###############################################################################
# STEP 3: do some post-processing
###############################################################################

# replace low-frequency codons
# -----------------------------------------------------------------------------
with open("results/optimized_sequence_for_poae.txt" , mode = 'r') as f:
    dnaSeq = f.read().strip()
cd = cas.CodingDomain(dnaSeq)   # read optimized sequence from text file
print(cd.as_DNAstring())

# decide which codons are too low in frequency and find replacement
replacements = cd.get_codonFrequency().lowFrequentCodonsAndReplacements(0.1, normalized = True)
# do the replacement
cd.replaceCodons(replacements)  

# replace unwanted sequencies (typically restriction sites)
# -----------------------------------------------------------------------------
extractor = cas.CodonFrequencyExtractor('data/highly_expressed_cds.txt')
   
# transform into AA sequence
AAseq = cas.translateDNA(dnaSeq)
optimizer = cas.CodonOptimizer(AAseq, extractor.codonFrequency ,extractor.codonPairFrequency)
optimizer.aaSequence.set_codingDomain(cd)


# list strings that cannot be present
restriction_sites = ["GGTACC", "TTCGAA", "TTAATTAA", "CACGTG", "AGATCT", "ACTAGT", "CTCGAG", "GTCGAC", "AAGCTT", "GAGCTC", "GGATCC"]

# let the optimizer find optimal replacements
optimizer.replaceUnwantedSubstrings(restriction_sites)
cd.get_codonFrequency()
# FINALLY: look at the result
optimizer.get_codingDomain().as_DNAstring()


    
    