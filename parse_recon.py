import scipy.io
import numpy as np
from itertools import combinations
import argparse

parser = argparse.ArgumentParser(description='Run Crossvalidation Pipeline.')

parser.add_argument('--input', "-i", type=str, default="Recon3D_301/Recon3D_301.mat",
                    help='Path to the Recon3D matlab file.')
parser.add_argument('--output', "-o", type=str, default="recon.tsv",
                    help='Path to store the edgelist file.')
parser.add_argument('--blacklist', "-b", type=str, default="reconparser/blacklist_metabolites.txt",
                    help='Path to store the edgelist file.')
parser.add_argument('--filtered', "-f", type=str, default="filtered_metabolites.txt",
                    help='Path to store filtered metabolties list.')



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
for i, col in enumerate(range(genes2metabolites.shape[1])):
    geneset = entrez_ids[genes2metabolites[:, col]].tolist()
    metabolite = metabolites[i][0][0]
    if len(geneset) > 0:
        edges = list(combinations(geneset, 2))
        written_edges += len(edges)
        with open(output_file, "a") as edgelist:
            for edge in edges:
                edgelist.write("{}\t{}\t{}\n".format(edge[0][0][0].split(".")[0], edge[1][0][0].split(".")[0], metabolite))

    if i % 500 == 0:
        print("Parsed {} out of {} metabolites. Wrote {} Edges.".format(i, genes2metabolites.shape[1], written_edges))