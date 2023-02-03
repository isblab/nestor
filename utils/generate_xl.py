import os, sys
import math
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdb_file", type=str, required=True, help="Input PDB file")
parser.add_argument("--dt", type=float, required=True, help="Linker max length")
parser.add_argument("--fdr", type=float, default=0.05, help="Linker max length")
parser.add_argument("--total_xl", type=int, required=True, help="No. of XL")
parser.add_argument("--split", type=float, default=0.7, help="Prior-evicalc split")
args = parser.parse_args()


class particle:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


def get_distance(particle1, particle2):
    raw_d = (
        (particle1.x - particle2.x) ** 2
        + (particle1.y - particle2.y) ** 2
        + (particle1.z - particle2.z) ** 2
    )
    d = math.sqrt(raw_d)
    return d


def select_random_entries(input_list, num_out):
    random.shuffle(input_list)
    outlist = []
    while len(outlist) < num_out:
        choice = random.choice(input_list)
        if choice not in outlist:
            outlist.append(choice)
    return outlist


def write_xl_file(fname, xl_list):
    with open(fname, "w") as outf:
        for lnk in xl_list:
            outf.write(lnk + "\n")


def create_pseudobond_file(fname, xl_list, color):
    model_num = 0
    atom_type = "ca"
    with open(fname, "w") as outf:
        for link in xl_list:
            chain1 = link.split(",")[0][-1].strip()
            res1 = int(link.split(",")[1])
            chain2 = link.split(",")[2][-1].strip()
            res2 = int(link.split(",")[3])
            pseudobond = f"#{model_num}:{res1}.{chain1}@{atom_type} #{model_num}:{res2}.{chain2}@{atom_type} {color}"
            outf.write(pseudobond + "\n")


##########################################################################################################################################
################################################################### Main #################################################################
##########################################################################################################################################

pdb_lines = {}
with open(args.pdb_file, "r") as pdb:
    for ln in pdb.readlines():
        ln = ln.strip()
        if ln.startswith("ATOM"):
            if ln[13:15] == "CA":
                coordinates = particle(
                    float(ln[31:38]), float(ln[39:46]), float(ln[46:55])
                )
                pdb_lines[coordinates] = (ln, len(pdb_lines.keys()))

contacts = []
non_contacts = []
for coords1, p1 in pdb_lines.items():
    for coords2, p2 in pdb_lines.items():
        if coords1 != coords2:
            if abs(p1[1] - p2[1]) > 3:
                distance = get_distance(coords1, coords2)
                if distance <= args.dt:
                    contacts.append((p1[0], p2[0]))
                else:
                    non_contacts.append((p1[0], p2[0]))

num_false_xl = int(args.fdr * args.total_xl)
num_true_xl = int(args.total_xl - num_false_xl)
false_xls = select_random_entries(non_contacts, num_false_xl)
xlinked_residues = select_random_entries(contacts, num_true_xl)
print(len(false_xls), len(xlinked_residues))

for lnk in false_xls:
    xlinked_residues.append(lnk)
random.shuffle(xlinked_residues)

prior_xlinks = []
li_xlinks = []
for entry in xlinked_residues:
    link = (
        f"chain{entry[0][21]},{entry[0][22:26]},chain{entry[1][21]},{entry[1][22:26]}"
    )

    if random.random() <= args.split:
        li_xlinks.append(link)
    else:
        prior_xlinks.append(link)


write_xl_file("prior_pseudo_xl.dat", prior_xlinks)
write_xl_file("evicalc_pseudo_xl.dat", li_xlinks)

create_pseudobond_file("prior_pseudobond.txt", prior_xlinks, "#00bfff")
create_pseudobond_file("evicalc_pseudobond.txt", li_xlinks, "#ea4335")
