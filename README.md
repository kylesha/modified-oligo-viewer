# modified-oligo-viewer
#### Python 3 class to facilitate visualization of oligonucleotides containing modified bases.

## Background
Molecular biologists often use modified bases in oligonucleotides (aka "oligos") to change the properties of the oligo in certain applications. For example, an oligo containing one or more "Super-T" bases has a higher melting temperature (Tm) compared to an oligo of identical sequence containing only canonical (natural) T bases (https://www.idtdna.com/site/Catalog/Modifications/Category/7). Thus, Super-T can be incorporated into probes or PCR primers that target AT-rich regions, to compensate for the lower Tm.

## The Problem
However, visualizing oligos containing modified bases can be a clumsy process. For example, suppose we have the oligo `ATCGAGTTTACCATATCTAGAATGCAT` with the 7th T replaced by the modified base "Super T" from [Integrated DNA Technologies (IDT)](https://www.idtdna.com/site/Catalog/Modifications/Category/7). A common visualization is to use a unique string representation of the modified base. For example,  would represent the modified oligo as `ATCGAG/iSuper-dT/TTACCATATCTAGAATGCAT`. Suppose we have another version of the oligo where the 7th and 9th T's are SuperT: `ATCGAG/iSuper-dT/T/iSuper-dT/ACCATATCTAGAATGCAT`. And another version where the first A is replaced by ["Fluoro A"](https://www.idtdna.com/site/Catalog/Modifications/Product/2791). Now we have a collection of slightly different versions of the same oligo:

`oligo1: ATCGAGTTTACCATATCTAGAATGCAT`
`oligo2: ATCGAG/iSuper-dT/TTACCATATCTAGAATGCAT`
`oligo3: ATCGAG/iSuper-dT/T/iSuper-dT/ACCATATCTAGAATGCAT`
`oligo4: /52FA/ATCGAGTTTACCATATCTAGAATGCAT`

If there many versions of the same oligo, it can be difficult to visualize and compare them on a base position basis because they are out of alignment. This clumsiness is compounded if the oligo has a longer sequence and contains multiple (different) modifications on each oligo.

## A Solution
One (of many) solution to better visualize a collection of modified oligos is to
1. allow the user to define a single-character symbol to represent the modified base. This resolves the issue of alignment.
2. colorize each of the four canonical bases (A, C, G, T). Additionally, assign the same color to the symbol representation of each base.


## Usage
The `ModifiedOligo` class was written to achieve these objectives.

__See the Jupyter notebook for an example use case.__

#### Dependencies
The `ModifiedOligo` class requires the [Biopython package](https://biopython.org/wiki/Download). It uses Biopython's `Bio.SeqRecord` object as the starting point for creating an oligo. The rational is that the `ModifiedOligo` object can then be added back as metadata to BioPython's `Bio.SeqFeature` object. (However, this feature is not yet implemented at this time). Additionally, Biopython enforces creation of unambiguous biological sequences. So for example, "AGTA" is a DNA sequence of adenine+guanine+thymine+adenine, NOT a peptide of alanine+glycine+threonine+alanine!  

After installation, perform the following imports:

```
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
```

#### Creating a `ModifiedOligo` object
To create a `ModifiedOligo` object, supply two objects to the class constructor:
- a `SeqRecord` object.
- a dictionary of the form {(symbol, name):(positions)}. The dictionary maps a user-defined symbol and the modified base name to the base position(s) in which the modified base occurs. The base numbering is __zero-based__ (so as to be consistent with `SeqRecord` base numbering convention).

For example, suppose we have an oligo (say `ATGTCAGTC`) in which the 2nd and 8th T are deoxyUracil (dU) bases, and the 9th C is a dideoxy C (ddC). We wish to represent dU using '#' and ddC using '$', respectively.

```
# construct IUPAC DNA object
dna = Seq("ATGTCAGTC", IUPAC.unambiguous_dna)

# construct SeqRecord object
seq_rec = SeqRecord(dna)  

# define modified bases, their symbol representations, and their positions in the oligo
mods = {('#', 'dU'):(1, 7),
        ('$', "ddC"):(8, )}   # singletone tuple MUST end with comma

# create ModifiedOligo object
mod_oligo = ModifiedOligo(seq_rec, modifications=mods)
```

#### Viewing a `ModifiedOligo` object
`ModifiedOligo` has two instance methods with the following signature:
- `view53(modified=False, showLegend=False)`: view the oligo in the `5'->3'` direction
- `view35(modified=False, showLegend=False)`: view the oligo in the `3'->5'` direction

Each method will pad the appropriate 5' or 3' to each end of the oligo. When `modified=False` (default), the unmodified (original) version of the oligo is displayed. To display the modified version, set `modified=False`. When the `showLegend` argument is True, it will show a legend indicating the symbol and the modified base that it represents.
