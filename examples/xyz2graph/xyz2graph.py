import re
import string
from itertools import combinations
from math import sqrt
from rdkit import Chem
from pathlib import Path
from rdkit.Chem import AllChem

atomic_radii = dict(Ac=1.88, Ag=1.59, Al=1.35, Am=1.51, As=1.21, Au=1.50, B=0.83, Ba=1.34, Be=0.35, Bi=1.54, Br=1.21,
                    C=0.68, Ca=0.99, Cd=1.69, Ce=1.83, Cl=0.99, Co=1.33, Cr=1.35, Cs=1.67, Cu=1.52, D=0.23, Dy=1.75,
                    Er=1.73, Eu=1.99, F=0.64, Fe=1.34, Ga=1.22, Gd=1.79, Ge=1.17, H=0.23, Hf=1.57, Hg=1.70, Ho=1.74,
                    I=1.40, In=1.63, Ir=1.32, K=1.33, La=1.87, Li=0.68, Lu=1.72, Mg=1.10, Mn=1.35, Mo=1.47, N=0.68,
                    Na=0.97, Nb=1.48, Nd=1.81, Ni=1.50, Np=1.55, O=0.68, Os=1.37, P=1.05, Pa=1.61, Pb=1.54, Pd=1.50,
                    Pm=1.80, Po=1.68, Pr=1.82, Pt=1.50, Pu=1.53, Ra=1.90, Rb=1.47, Re=1.35, Rh=1.45, Ru=1.40, S=1.02,
                    Sb=1.46, Sc=1.44, Se=1.22, Si=1.20, Sm=1.80, Sn=1.46, Sr=1.12, Ta=1.43, Tb=1.76, Tc=1.35, Te=1.47,
                    Th=1.79, Ti=1.47, Tl=1.55, Tm=1.72, U=1.58, V=1.33, W=1.37, Y=1.78, Yb=1.94, Zn=1.45, Zr=1.56)


class MolGraph:
    """Represents a molecular graph."""

    def __init__(self):
        self.elements = []
        self.x = []
        self.y = []
        self.z = []
        self.adj_list = {}
        self.atomic_radii = []
        self.bond_lengths = {}
        self.mapping_list = {}
        self.offset = []
        self.initial_state = []
        self.end_state = []
        self.atom_idx = []

    def read_xyz(self, file_path: str) -> None:
        """Reads an XYZ file, searches for elements and their cartesian coordinates
        and adds them to corresponding arrays."""
        pattern = re.compile(r'([A-Za-z]{1,3})\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)\s*(-?\d+(?:\.\d+)?)')
        with open(file_path) as file:
            for element, x, y, z in pattern.findall(file.read()):
                self.elements.append(element)
                self.x.append(float(x))
                self.y.append(float(y))
                self.z.append(float(z))
        self.atomic_radii = [atomic_radii[element] for element in self.elements]
        self._generate_adjacency_list()
    
    def draw_mol(self, mol, pos_offset=0):
        pos = mol.GetConformer().GetPositions() 
        pos[:, 0] = pos[:, 0] - pos_offset

        for i in range(mol.GetNumAtoms()):
            self.elements.append(mol.GetAtomWithIdx(i).GetSymbol())
        self.x += pos[:, 0].tolist()
        self.y += pos[:, 1].tolist()
        self.z += pos[:, 2].tolist()

        self.atomic_radii = [atomic_radii[element] for element in self.elements]
        self.generate_adj(mol)
        self.offset.append(mol.GetNumAtoms())

    def read_smiles(self, smiles, pos_offset=0):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        self.draw_mol(mol, pos_offset)

    def generate_mapping(self, mcs):
        for pair in mcs:
            i, j = pair
            j = j + self.offset[-2]
            self.mapping_list.setdefault(i, set()).add(j)
            self.mapping_list.setdefault(j, set()).add(i)

    def read_molfile(self, molfile, pos_offset=0):
        mol = Chem.MolFromMolFile(molfile, removeHs=False)
        self.draw_mol(mol, pos_offset)

    def read_grofile(self, gro_file):
        content = Path(gro_file).read_text()
        top_file = Path(gro_file).with_suffix('.top')
        self.initial_state, self.end_state = self.parse_topfile(top_file)

        count = 0
        for row in content.split('\n'):
            if '1LIG' in row:
                new_row = self.strip_blank(row)
                self.atom_idx.append(new_row[2])
                pos = list(map(float, new_row[3:]))
                self.elements.append(self.initial_state[count])
                self.x.append(pos[0]*10)
                self.y.append(pos[1]*10)
                self.z.append(pos[2]*10)
                count += 1
        for idx, element in enumerate(self.elements):
            if element == '*':
                to_ele = self.end_state[idx]
                self.atomic_radii.append(atomic_radii[to_ele])
            else:
                self.atomic_radii.append(atomic_radii[element])
        self._generate_adjacency_list()

    def strip_blank(self, row):
        new_row = []
        for item in row.split(' '):
            if item:
                new_row.append(item)
        return new_row

    def parse_topfile(self, top_file):
        content = Path(top_file).read_text()
        peroidic_table = Chem.GetPeriodicTable()
        name_to_symbol = {}
        atom_type_field = False
        atom_pair_field = False
        initial_state = []
        end_state = []
        for row in content.split('\n'):
            if 'at.num' in row:
                atom_type_field = True
                continue
            if atom_type_field:
                if not row:
                    atom_type_field = False
                    continue
                new_row = self.strip_blank(row)
                name_to_symbol[new_row[0]] = peroidic_table.GetElementSymbol(int(new_row[1]))
            if 'type1' in row:
                atom_pair_field = True
                continue
            if atom_pair_field:
                if not row:
                    atom_pair_field = False
                    continue
                new_row = self.strip_blank(row)
                initial_state.append(name_to_symbol[new_row[1]])
                end_state.append(name_to_symbol[new_row[8]])
        return initial_state, end_state

    def generate_adj(self, mol):
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if self.offset:
                i, j = i + self.offset[-1], j + self.offset[-1]

            x_i, y_i, z_i = self.__getitem__(i)[1]
            x_j, y_j, z_j = self.__getitem__(j)[1]
            distance = sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)

            self.adj_list.setdefault(i, set()).add(j)
            self.adj_list.setdefault(j, set()).add(i)
            self.bond_lengths[frozenset([i, j])] = round(distance, 5)

    def _generate_adjacency_list(self):
        """Generates an adjacency list from atomic cartesian coordinates."""
        node_ids = range(len(self.elements))
        for i, j in combinations(node_ids, 2):
            x_i, y_i, z_i = self.__getitem__(i)[1]
            x_j, y_j, z_j = self.__getitem__(j)[1]
            distance = sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)
            if 0.1 < distance < (self.atomic_radii[i] + self.atomic_radii[j]) * 1.3:
                self.adj_list.setdefault(i, set()).add(j)
                self.adj_list.setdefault(j, set()).add(i)
                self.bond_lengths[frozenset([i, j])] = round(distance, 5)


    def edges(self):
        """Creates an iterator with all graph edges."""
        edges = set()
        for node, neighbours in self.adj_list.items():
            for neighbour in neighbours:
                edge = frozenset([node, neighbour])
                if edge in edges:
                    continue
                edges.add(edge)
                yield node, neighbour

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, position):
        return self.elements[position], (
            self.x[position], self.y[position], self.z[position])
