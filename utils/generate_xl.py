import os, sys
import math
import random

pdb_file = sys.argv[1]
distance_threshold = float(sys.argv[2])
total_xl = int(sys.argv[3])
split = float(sys.argv[4])


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


def create_pseudobond_file(fname, xl_list):
    model_num = 0
    atom_type = "ca"
    with open(fname, "w") as outf:
        for link in xl_list:
            chain1 = link.split(",")[0][-1]
            chain2 = link.split(",")[2][-1]
            pseudobond = f"#{model_num}:{link.split(',')[1]}.{chain1}@{atom_type} #{model_num}:{link.split(',')[3]}.{chain2}@{atom_type} #00bfff"
            outf.write(pseudobond + "\n")


##########################################################################################################################################
################################################################### Main #################################################################
##########################################################################################################################################

ca_particles = []
pdb_lines = []
with open(pdb_file, "r") as pdb:
    for ln in pdb.readlines():
        ln = ln.strip()
        if ln.startswith("ATOM"):
            if ln[13:15] == "CA":
                ca_particles.append(particle(ln[31:38], ln[39:46], ln[46:55]))
                pdb_lines.append(ln)

contacts = []
for idx1, p1 in enumerate(ca_particles):
    for idx2, p2 in enumerate(ca_particles):
        if idx1 != idx2:
            if abs(idx1 - idx2) > 3:
                distance = get_distance(p1, p2)
                if distance <= distance_threshold:
                    contacts.append((idx1, idx2))
                    # print(f"{pdb_lines[idx1]}{pdb_lines[idx2]}\n")


xlinked_atom_ids = select_random_entries(contacts, total_xl)
prior_xlinks = []
li_xlinks = []

for entry in xlinked_atom_ids:
    linked_residues = (pdb_lines[entry[0]], pdb_lines[entry[1]])
    link = f"chainA,{linked_residues[0][22:26]},chainA,{linked_residues[1][22:26]}"

    if random.random() <= split:
        li_xlinks.append(link)
    else:
        prior_xlinks.append(link)

write_xl_file("prior_pseudo_xl.dat", prior_xlinks)
write_xl_file("evicalc_pseudo_xl.dat", li_xlinks)

create_pseudobond_file("prior_pseudobond.txt", prior_xlinks)
create_pseudobond_file("evicalc_pseudobond.txt", li_xlinks)
