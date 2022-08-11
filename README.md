# reconparser
A Parser for Recon3D that produces Gene Edgelists

This parser produces edgelists where two genes are connected if they share a metabolite. Future versions will include directionality and/or weights.

To avoid very large clusters, it is possible and recommended to provide a blacklist to filter out very common metabolites like h2o.
The blacklist is used as prefixes, so ```h2o``` would also filter ```h2o2```, if you want to prevent that explicitely enter ```h20[``` with the opening square bracket.

## Download Recon3D

```
./download.sh
```

## Run Parser

```
python3 parse_recon -i <path to recon .mat file> -o <path and filename of output edgelist> -b [OPTIONAL]<path to blacklist> -f [OPTIONAL]<path to list of filtered metabolties>
```

For example just run it as follows:

```
python3 parse_recon -i Recon3D_301/Recon3D_301.mat -o recon.tsv -b blacklist_metabolites.tsv
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
