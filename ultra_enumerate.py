from rdkit import Chem
from itertools import permutations
import random
import re
from math import factorial
import argparse
import csv




def get_tokens(mol):
    tokens = []
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        token = sym
        if atom.GetIsAromatic():
            token = token.lower()
        hs = atom.GetTotalNumHs()
        iso = atom.GetIsotope()
        cha = atom.GetFormalCharge()
        if cha or iso or token=="n":
            if hs:
                token += "H"
                if hs>1:
                    token += str(hs)
        if iso:
            token = str(iso)+token
        if cha:
            if cha>0:
                token += "+"
                if cha>1:
                    token += str(cha)
            if cha<0:
                token += "-"
                if cha>1:
                    token += str(-cha)
        if cha or iso or sym not in ["B","C","N","O","F","Cl","Br","I"]:
            token = f'[{token}]'
        elif token=="nH":
            token = f'[{token}]'
        tokens.append(token)
    return tokens
    
float2bond = {1:"",1.5:":",2:"=",3:"#"}

def get_bonds(m):
    bondtypes = {}
    bonds = []
    am = Chem.GetAdjacencyMatrix(m,useBO=True)
    for i,row in enumerate(am):
        for j,el in enumerate(row):
            if i<j:
                if am[i][j]:
                    bondtypes[(i,j)] = float2bond[am[i][j]]
                    bonds.append((i,j))
    return bonds,bondtypes

def get_bondtokens(bonds,bondtypes,num_atoms):
    bondtokens = [""]*num_atoms
    curr_idx = 1
    for i,bond in enumerate(bonds):
        bondtokens[bond[0]]+=bondtypes[bond]
        bondtokens[bond[1]]+=bondtypes[bond]
        bondtokens[bond[0]]+="("
        bondtokens[bond[1]]+="("
        bondtokens[bond[0]]+= str(curr_idx)
        bondtokens[bond[1]]+= str(curr_idx)
        bondtokens[bond[0]]+=")"
        bondtokens[bond[1]]+=")"
        curr_idx += 1
    return bondtokens

def remap_mol(mappings, mol):
    atom_tokens = get_tokens(mol)
    bonds, bondtypes = get_bonds(mol)
    bond_tokens = get_bondtokens(bonds,bondtypes,len(atom_tokens))
    smis = []
    for mapping in mappings:
        all_tokens = [atom_tokens[idx]+bond_tokens[idx] for idx in mapping]
        smi = ".".join(all_tokens)
        matches = re.findall(r"\([^()]*\)", smi)
        final_tokens = {}
        curr_idx = 1
        for match in matches:
            if match not in final_tokens:
                if curr_idx>9:
                    final_tokens[match] = f"%{curr_idx}"
                else:
                    final_tokens[match] = f"{curr_idx}"
                curr_idx+=1
        for k in final_tokens:
            smi = smi.replace(k,final_tokens[k])
        smis.append(smi)
    return smis
    
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Generate unique random permutations of [0..N] with SMILES input."
    )
    
    parser.add_argument(
        "-N", 
        type=int, 
        default=100,
        help="Size of the permutation set (default: 100)"
    )
    
    parser.add_argument(
        "--smiles",
        type=str,
        default="c1ccccc1",
        help="SMILES string (default: c1ccccc1)"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output file (default: None, prints to stdout)"
    )
    
    args = parser.parse_args()
    
    smi = args.smiles
    m = Chem.MolFromSmiles(smi)
    smis = []
    atomcount = len(m.GetAtoms())
    fac = factorial(atomcount)
    perms = set([])
    N = args.N
    N = min(N,fac)
    while len(perms)<N:
        perms.add(tuple(random.sample(list(range(atomcount)),atomcount)))
    smis = remap_mol(perms,m)
    if args.output:
        with open(args.output, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=' ',
                                    quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for smi in smis:
                writer.writerow([smi])


    else:
        for smi in smis:
            print(smi)

if __name__ == "__main__":
    main()

