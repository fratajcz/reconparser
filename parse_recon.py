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
parser.add_argument('--use_direction', "-d", action='store_true',
                    help="Adds an edge from A to B only if A produces what B consumes. NOT IMPLEMENTED")
parser.add_argument('--use_weight', "-w", action='store_true',
                    help="Adds an additional column 'Weight' which displays stochiometry. NOT IMPLEMENTED")

args = parser.parse_args()


def generate_edges(producing_geneset, consuming_geneset=None):
    """ generates edges from genesets. 
    
    If only producing geneset is provided, it is assumed that we are in undirected setting and all possible combinations are produces
    If both sets are provided, only edges running from all producing to all consuming genes are provided.

    RETURNS list of tuples
    
    """

    if consuming_geneset is None:
        edges = list(combinations(geneset, 2))
    else:
        edges = [[producer, consumer] for consumer in consuming_geneset for producer in producing_geneset]
    return edges


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

if not args.use_direction:
    genes2metabolites[genes2metabolites != 0] = 1
    genes2metabolites = genes2metabolites.astype(np.bool8)
else:
    producing_genes2metabolites = genes2metabolites.copy().astype(np.bool8)
    consuming_genes2metabolites = genes2metabolites.copy().astype(np.bool8)

    producing_genes2metabolites[genes2metabolites > 0] = 1
    producing_genes2metabolites[genes2metabolites <= 0] = 0

    consuming_genes2metabolites[genes2metabolites < 0] = 1
    consuming_genes2metabolites[genes2metabolites >= 0] = 0



entrez_ids = recon[9]

with open(output_file, "w") as edgelist:
    edgelist.write("{}\t{}\t{}\n".format("EntrezA", "EntrezB", "Metabolite"))

written_edges = 0

if args.ignore_components:
    all_edges = set()

for i, col in enumerate(range(genes2metabolites.shape[1])):

    if not args.use_direction:
        geneset = entrez_ids[genes2metabolites[:, col]].tolist()
        not_empty = len(geneset) > 0
        genesets = [geneset]
    else:
        producing_geneset = entrez_ids[producing_genes2metabolites[:, col]].tolist()
        consuming_geneset = entrez_ids[consuming_genes2metabolites[:, col]].tolist()
        genesets = [producing_geneset, consuming_geneset]

        not_empty = (len(producing_geneset)) > 0 and (len(consuming_geneset) > 0)

    metabolite = metabolites[i][0][0]
    metabolite = metabolite.split("[")[0] if args.ignore_components else metabolite

    if not_empty:
        edges = generate_edges(*genesets)
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
