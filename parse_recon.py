import scipy.io
import numpy as np
from itertools import combinations
import argparse

parser = argparse.ArgumentParser(description='Run Crossvalidation Pipeline.')

parser.add_argument('--input', "-i", type=str, default="../Recon3D_301/Recon3D_301.mat",
                    help='Path to the Recon3D matlab file.')
parser.add_argument('--output', "-o", type=str, default="data/recon.tsv",
                    help='Path to store the edgelist file.')
parser.add_argument('--blacklist', "-b", type=str, default="blacklist_metabolites.txt",
                    help='Path to store the edgelist file.')
parser.add_argument('--filtered', "-f", type=str, default="data/filtered_metabolites.txt",
                    help='Path to store filtered metabolties list.')
parser.add_argument('--ignore-components', "-c", action='store_true',
                    help="Wheter component codes ([c], [l], etc.) should be ignored.")



args = parser.parse_args()


mat = scipy.io.loadmat(args.input)
recon = mat["Recon3D"][0][0]

output_file = args.output
filtered_list = args.filtered

metabolite_blacklist = []

with open(args.blacklist, "r") as file:
    for line in file:
        metabolite_blacklist.append(line.strip())

print("Read {} banned metabolite prefixes from {}.".format(len(metabolite_blacklist), args.blacklist))

metabolites = recon[1]

mask = np.ones_like(metabolites.flatten(), dtype=np.bool8)

for i, metabolite in enumerate(metabolites):
    for banned_metabolite in metabolite_blacklist:
        if metabolite[0][0].startswith(banned_metabolite):
            mask[i] = False

metabolites2reactions = recon[0][mask, :]
filtered_metabolites = metabolites[~mask]
metabolites = metabolites[mask]
reactions2genes = recon[21]

print("filtered out {} from a total of {} metabolites. See {} for a list of them.".format(np.sum(~mask), len(mask), filtered_list))

with open(filtered_list, "w") as file:
    for metabolite in filtered_metabolites:
        file.writelines("{}\n".format(metabolite[0][0]))

genes2metabolites = (metabolites2reactions * reactions2genes).transpose()

genes2metabolites[genes2metabolites != 0] = 1
genes2metabolites = genes2metabolites.astype(np.bool8)

entrez_ids = recon[9]

with open(output_file, "w") as edgelist:
    edgelist.write("{}\t{}\t{}\n".format("EntrezA", "EntrezB", "Metabolite"))

written_edges = 0

if args.ignore_components:
    all_edges = set()
for i, col in enumerate(range(genes2metabolites.shape[1])):
    geneset = entrez_ids[genes2metabolites[:, col]].tolist()
    metabolite = metabolites[i][0][0]
    metabolite = metabolite.split("[")[0]
    if len(geneset) > 0:
        edges = list(combinations(geneset, 2))
        written_edges += len(edges)
        with open(output_file, "a") as edgelist:
            for edge in edges:
                edge_tuple = (edge[0][0][0].split(".")[0], edge[1][0][0].split(".")[0], metabolite)
                if args.ignore_components:
                    all_edges.add(edge_tuple)
                else:
                    edgelist.write("{}\t{}\t{}\n".format(*edge_tuple))

    if i % 500 == 0 or i == genes2metabolites.shape[1] - 1:
        verb = "Gathered" if args.ignore_components else "Wrote"
        print("Parsed {} out of {} metabolites. {} {} Edges.".format(i, genes2metabolites.shape[1], verb, written_edges))

if args.ignore_components:
    with open(output_file, "a") as edgelist:
        all_edges = list(all_edges)
        all_edges.sort(key=lambda x: x[2])
        
        for i, edge in enumerate(all_edges):
            edgelist.write("{}\t{}\t{}\n".format(*edge))

            if i % 100000 == 0 or i == len(all_edges) - 1:
                print("Written {} out of {} deduplicated Edges.".format(i, len(all_edges)))
