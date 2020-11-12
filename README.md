# cas9optimization

`cas9optimization` is Python package (single module, pure Python) that can be used perform codon optimization. The coding domain of a protein is optimized such that the frequency distribution of the codon pairs matches that of the distribution that is observed in highly expressed coding domains of the target organism.

Even though the package is generic, it was designed to codon-optimize the CAS9 coding domain for the paper [1]

### Methodological aspects

The methodology is mainly borrowed from [2], where the coding domain of a protein is host-optimized by: 

1. **estimating** a set of **statistics** on the preferred **codon usage** of the host organism based on highly expressed protein-coding domains
2. **optimize** the **coding domain** of a target protein such that the codon statistics match those estimated above.
3. we added post **post-processing** functions (not described in [2]),  to avoid restriction sites and replace codons with a very low frequency

The implementation allows to use two sets of statistics: (a) the (traditional) observed codon frequency distribution; and (b) the observed codon-pair frequency distribution, where a codon-pair consists of two consecutive codons.

The main component of this package is an optimized *Simulated Annealing* (SA) algorithm that optimizes the coding domain (whereas [2] uses a genetic algorithm). Quite some work was invested in reducing the runtime of the algorithm by applying some computational tricks that allow to recompute the similarity between the host's codon-pair frequency distribution and that of the target protein after making a small change to the coding domain of the target protein. As a result the runtime is about 5 sec for a target protein such as CAS9. Note that his algorithm remains heuristic (and stochastic) and does not guarantee globally optimal solutions. 

### Examples

An extensive example can be found in `cas9_optimization_main.py` . Below, you can find the key-steps.

As a **first step**, a file containing the several highly expressed coding domains, one coding domain (DNA) string per line. In the `./data` folder, the file `highly_expressed_cds.txt` contains a set of highly expressed coding domains for *Fusarium poae*.  Lines starting with \# are ignored. These files are loaded and processed by a `CodonFrequencyExtractor`.

 ```python
import cas9_optimization as cas
extractor = cas.CodonFrequencyExtractor('data/highly_expressed_cds.txt')
print(cas.AA_TO_CODON)           # dictionary of AAs (keys) and related codons (values)
print(extractor.codonFrequency)  # observed codon frequencies
print(extractor.codonPairFrequency)  # observed codon pair frequencies
 ```

(Optional) The frequencies can be interpreted by combining `cas.AA_TO_CODON` and the observed frequencies.

```python
>>> print(cas.AA_TO_CODON['G'])     # codons coding for AA with single letter code G
['GGT', 'GGC', 'GGA', 'GGG']
>>> extractor.codonFrequency['G']   # observed frequencies for these codons
array([536, 288,  92,  14])
```

```python
>>> print(cas.AAPAIR_TO_CODON['GF'])     # codon pairs encoding for AA pair GF (SLCs)
['GGTTTT', 'GGTTTC', 'GGCTTT', 'GGCTTC', 'GGATTT', 'GGATTC', 'GGGTTT', 'GGGTTC']
>>> extractor.codonPairFrequency['GF']
array([ 0, 23,  3,  4,  1,  0,  0,  0])
```

As a **second step**, the AA sequence of the protein that should be imported as a string. The file `cas9_coding_domain_T_reesei.txt` contains the CAS9 optizimized for *T. reesei* that is imported and translated into a DNA sequence.

```python
with open("data/cas9_coding_domain_T_reesei.txt", mode = "r") as f:
    DNAseq = f.readlines()[1].strip()
# transform into AA sequence
AAseq = cas.translateDNA(DNAseq)
```

As a **third step**, an `CodonOptimizer` object is created and applied.

```python
# create an optimizer object
optimizer = cas.CodonOptimizer(AAseq, extractor.codonFrequency,
                                   extractor.codonPairFrequency)
# perform codon optimization where the frequency distribution of codon pairs
# is choosen to match the frequencey observed in the highly expressed domains
message = optimizer.optimizeCodingDomainBivariate(algorithm="SA")

# extract codon-optimzed DNA string
print(optimizer.get_codingDomain().as_DNAstring())
```

As a **final step** ***postprocessing*** and ***quality checks*** can be applied (see `cas9_optimization_main.py` for examples )

### How to cite

If you use this code, please cite ... *(to be added)*

### Literature

[1] ... *(ours, to be added)*

[2]