# selectClassSite

`selectClassSite` is program to extract sites from a alignment of coding sequence in fasta or phylip format.
This program is written using the [**Bio ++** library ](https://biopp.github.io/) (Guéguen et al. 2013).

You can use the static executable compiled for linux x64 computer (see [Release](https://github.com/benoitnabholz/selectClassSite/releases/)). You can also compile the program assuming that you have [**Bio ++**](https://biopp.github.io/) installed (here the Bio++ library V2 is locally installed in `$HOME/local/bpp/dev/` directory).

Please cite Bio++ if you use this program:
- Guéguen L, Gaillard S, Boussau B, Gouy M, Groussin M, Rochette NC, Bigot T, Fournier D, Pouyet F, Cahais V, et al. 2013. Bio++: Efficient Extensible Libraries and Tools for Computational Molecular Evolution. Mol. Biol. Evol. 30:1745–1750.

*Author:* Benoit Nabholz

## Option

selectClassSite file_name fomat 'Standard or VertebrateMitochondrial or InvertebrateMitochondrial' 'L4 or third or independent_codon' 

- format : fasta or phylip
- third: Extract third codon position
- L4: Extract fourfold degenerate site
- independent_codon: print two files with independent codon (selected randomly)

### Example of usage to extract the 4-fold denerate site from a alignment (align.fata)
`selectClassSite align.fasta fasta Standard L4`
