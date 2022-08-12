# reconparser
A Parser for Recon3D that produces Gene Edgelists

This parser produces edgelists where two genes are connected if they share a metabolite. Future versions will include directionality and/or weights.

To avoid very large clusters, it is possible and recommended to provide a blacklist to filter out very common metabolites like ```h2o```.
The blacklist is used as prefixes, so ```h2o``` would also filter ```h2o2```, if you want to prevent that explicitely enter ```h2o[``` with the opening square bracket.

## Shortcut

If you are fine with the selection of blacklisted metabolites in ```blacklist_metabolites.tsv```, you can go ahead and directly use the edgelist that is stored at ```data/recon.tsv```.

## Dependencies

```
numpy
scipy
```

## 1. Download Recon3D

```
./download.sh
```

## 2. Run Parser

```
python3 parse_recon.py -i <path to recon .mat file> -o <path and filename of output edgelist> -b [OPTIONAL]<path to blacklist> -f [OPTIONAL]<path to list of filtered metabolties>
```

For example just run it as follows:

```
python3 parse_recon.py -i Recon3D_301/Recon3D_301.mat -o recon.tsv -b blacklist_metabolites.tsv
```

and enjoy your edgelist:

```
EntrezA	EntrezB	Metabolite
10840	2356	10fthf[c]
10840	160428	10fthf[c]
10840	4522	10fthf[c]
10840	2618	10fthf[c]
10840	2618	10fthf[c]
10840	100287639	10fthf[c]
10840	441024	10fthf[c]
...
```

## 3. Additional Options

### Ignore Cellular Component

If the ```-c/--ignore-component``` flag is passed, then the edges will be filtered for multiple occurrences of (GeneA, GeneB, Metabolite) if two such edges have the same Metabolite, but in different cellular components (such as ```h2o[c]``` and ```h2o[l]```). This slightly reduces the number of edges.

### Directed Network

If the ```-d/--use-direction``` flag is passed, then the network becomes directed with edges running only from GeneA to GeneB if GeneA produces a metabolite that GeneB consumes. This breaks up the large cliques that the undirected version produces and heavily reduces the number of edges.

## 4. Future Work

### Weighted Network

Weights of edges in the network could reflect the number of reactions in which a gene produces/consumes a metabolite, or the average stochiometry of a metabolite in these reactions or the balanced-ness of the interaction of two genes (i.e. weight is 1 if geneA produces 2 mols of metaboliteX per reaction and geneB consumes 2 mols of metaboliteX per reaction). Due to this ambiguity of meaning, the implementation of weighted edges is postponed.
