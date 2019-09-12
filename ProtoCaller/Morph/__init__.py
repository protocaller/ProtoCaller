import copy as _copy
import os as _os
import re as _re

import parmed as _pmd

# TODO: make Morph viable
# TODO: convert print statements to logs
__all__ = []


class Morph():
    def __init__(self, system1, system2):
        for res1, res2 in zip(system1.residues, system2.residues):
            if res1.name != res2.name:
                raise ValueError("Need the same type and sequence of residues for both systems.")
        self.system1 = _copy.deepcopy(system1)
        self.system2 = _copy.deepcopy(system2)
        self.merge()

    def merge(self):
        def dummify(atom):
            atom.atom_type.epsilon = 0
            atom.atom_type.charge = 0
            atom.atom_type.rmin = 0
            atom.atom_type.name = "du"
            atom.atom_type.atomic_number = 0
            atom.charge = 0
            atom.name = "du"
            atom.type = "du"
            return atom

        self.map = []

        N = 0
        for r, (res1, res2) in enumerate(zip(self.system1.residues, self.system2.residues)):
            incr = 0
            N += len(res1.atoms)

            map = {}

            for i, atom1 in enumerate(res1.atoms):
                for j, atom2 in enumerate(res2.atoms):
                    if atom1.xx == atom2.xx and atom1.xy == atom2.xy and atom1.xz == atom2.xz:
                        map[i] = j
                        break
                    if j == len(self.system2.residues[r]) - 1:
                        du = _copy.deepcopy(atom1)
                        dummify(du)
                        self.system2.add_atom(du, res2.name, res2.number)
                        map[i] = N + incr
                        incr += 1
                        break

            rev_map = lambda x: list(map.keys())[list(map.values()).index(x)]

            for j, atom2 in enumerate(res2.atoms):
                try:
                    rev_map(j)
                except ValueError:
                    # calling the property enables index update - this seemingly pointless line is needed
                    du = _copy.deepcopy(atom2)
                    dummify(du)
                    self.system1.add_atom(du, res1.name, res1.number)
                    self.system1.atoms[-1].idx
                    map[len(res1.atoms) - 1] = j

            for j, atom2 in enumerate(res2.atoms):
                # calling the property enables index update - these two seemingly pointless lines are needed
                atom2.idx
                atom2._idx = rev_map(j)
                atom2.idx

            res2.sort()

        self.system1.atoms.sort()
        self.system2.atoms.sort()

        # bonds - needed since parmed doesn't write nonbonded atoms
        bonds_1 = {tuple(sorted([x.atom1.idx, x.atom2.idx])): x.type for x in self.system1.bonds}
        bonds_2 = {tuple(sorted([x.atom1.idx, x.atom2.idx])): x.type for x in self.system2.bonds}

        common_keys = set(bonds_1.keys()).intersection(bonds_2.keys())
        unique_bonds_1 = set(bonds_1.keys()) - common_keys
        unique_bonds_2 = set(bonds_2.keys()) - common_keys

        for i, j in ((1, 2), (2, 1)):
            exec("bonds_i, bonds_j = bonds_%d, bonds_%d" % (i, j), locals(), globals())
            exec("unique_bonds_i = unique_bonds_%d" % i, locals(), globals())
            exec("system_j = self.system%d" % j, locals(), globals())
            for key in unique_bonds_i:
                val = bonds_i[key]
                if val not in system_j.bond_types:
                    bond_type = _pmd.BondType(val.k, val.req, system_j.bond_types)
                    system_j.bond_types.append(bond_type)
                    bonds_j[key] = system_j.bond_types[-1]
                else:
                    bonds_j[key] = system_j.bond_types[system_j.bond_types.index(val)]
            system_j.bonds = bonds_j
            system_j.bonds = _pmd.TrackedList([_pmd.Bond(system_j.atoms[m], system_j.atoms[n], type=t)
                                               for (m, n), t in sorted(bonds_j.items())])

        '''
        # angles - not necessarily needed
        angles_1 = {(min(x.atom1.idx, x.atom3.idx), x.atom2.idx, max(x.atom1.idx, x.atom3.idx)) : x.type
                    for x in self.system1.angles}
        angles_2 = {(min(x.atom1.idx, x.atom3.idx), x.atom2.idx, max(x.atom1.idx, x.atom3.idx)): x.type
                    for x in self.system2.angles}

        common_keys = set(angles_1.keys()).intersection(angles_2.keys())
        unique_angles_1 = set(angles_1.keys()) - common_keys
        unique_angles_2 = set(angles_2.keys()) - common_keys

        for i, j in ((1, 2), (2, 1)):
            exec("angles_i, angles_j = angles_%d, angles_%d" % (i, j), locals(), globals())
            exec("unique_angles_i = unique_angles_%d" % i, locals(), globals())
            exec("system_j = self.system%d" % j, locals(), globals())
            for key in unique_angles_i:
                val = angles_i[key]
                if val not in system_j.angle_types:
                    angle_type = _pmd.AngleType(val.k, val.theteq, system_j.angle_types)
                    system_j.angle_types.append(angle_type)
                    angles_j[key] = system_j.angle_types[-1]
                else:
                    angles_j[key] = system_j.angle_types[system_j.angle_types.index(val)]
            system_j.angles = angles_j
            system_j.angles = _pmd.TrackedList([_pmd.Angle(system_j.atoms[m], system_j.atoms[n], system_j.atoms[o], type=t)
                                              for (m, n, o), t in sorted(angles_j.items())])


        # dihedrals - not necessarily needed; not completely functional for double dihedrals
        dihedrals_1 = {(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx) : x.type for x in self.system1.dihedrals}
        dihedrals_2 = {(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx) : x.type for x in self.system2.dihedrals}

        common_keys = set(dihedrals_1.keys()).intersection(dihedrals_2.keys())
        unique_dihedrals_1 = set(dihedrals_1.keys()) - common_keys
        unique_dihedrals_2 = set(dihedrals_2.keys()) - common_keys

        for i, j in ((1, 2), (2, 1)):
            exec("dihedrals_i, dihedrals_j = dihedrals_%d, dihedrals_%d" % (i, j), locals(), globals())
            exec("unique_dihedrals_i = unique_dihedrals_%d" % i, locals(), globals())
            exec("system_j = self.system%d" % j, locals(), globals())
            added_zero_dihedral = False
            for key in unique_dihedrals_i:
                val = dihedrals_i[key]
                #single dihedral
                if isinstance(val, _pmd.DihedralType):
                    if not added_zero_dihedral:
                        zero_dihedral = _pmd.DihedralType(0, 0, 0, list=system_j.dihedral_types)
                        system_j.dihedral_types.append(zero_dihedral)
                        added_zero_dihedral = True
                    dihedrals_j[key] = zero_dihedral
                #several dihedrals
                elif isinstance(val, _pmd.DihedralTypeList):
                    dihedral_list = _pmd.DihedralTypeList(list=system_j.dihedral_types)
                    dihedral_list_temp = [_pmd.DihedralType(0, k, 0, list=dihedral_list) for k, x in enumerate(val)]
                    for elem in dihedral_list_temp:
                        dihedral_list.append(elem)
                    dihedrals_j[key] = dihedral_list

            system_j.dihedrals = dihedrals_j
            system_j.dihedrals = _pmd.TrackedList([_pmd.Dihedral(system_j.atoms[m], system_j.atoms[n],
                                                               system_j.atoms[o], system_j.atoms[p], type=t)
                                                  for (m, n, o, p), t in sorted(dihedrals_j.items())])
        '''

    def write(self, filename, *args, **kwargs):
        ext = filename.split(".")[-1]
        if ext == "gro":
            self._writeToGro(filename)
        elif ext == "top":
            self._writeToTop(filename, *args, **kwargs)

    def _writeToGro(self, filename):
        try:
            _os.remove(filename)
        except:
            pass

        self.system1.save(filename)

    def _writeToTop(self, filename, intermediate_files=False):
        filebase, ext = _os.path.splitext(filename)[0], filename.split(".")[-1]
        filename1 = "%s_1.%s" % (filebase, ext)
        filename2 = "%s_2.%s" % (filebase, ext)

        for name in (filename, filename1, filename2):
            try:
                _os.remove(name)
            except:
                pass

        self.system1.save(filename1)
        self.system2.save(filename2)

        top1 = open(filename1).readlines()
        top2 = open(filename2).readlines()

        if not intermediate_files:
            _os.remove(filename1)
            _os.remove(filename2)

        # editing the two tops into a shape convenient for manipulation
        for j in (1, 2):
            exec("top = top%d" % j, locals(), globals())
            top_temp = []
            section_temp = []
            curr_section = None
            for i, line in enumerate(top):
                line = line.strip().split(";")[0]
                if len(line) or i == len(top) - 1:
                    match = _re.search(r"\[ (\w+) \]", line)
                    if match is not None and match.group(1) != curr_section:
                        if curr_section is not None:
                            top_temp += [(curr_section, section_temp)]
                            section_temp = []
                        curr_section = match.group(1)
                        continue
                    section_temp += [line.split()]
                    if i == len(top) - 1:
                        top_temp += [(curr_section, section_temp)]
            exec("top_edit%d = top_temp" % j, locals(), globals())

        # merging the two tops
        top12 = []
        for (sec1, list1), (sec2, list2) in zip(top_edit1, top_edit2):
            assert sec1 == sec2
            if list1 == list2:
                top12 += [(sec1, list1)]
                continue

            if sec1 == "defaults":
                print("Different defaults. Using the defaults from molecule 1...")
                top12 += [(sec1, list1)]
            elif sec1 == "atomtypes":
                atomtypes = {ls[0]: ls[1:] for ls in list1 + list2}
                top12 += [(sec1, [[key] + value for key, value in atomtypes.items()])]
            elif sec1 == "moleculetype":
                print("Different molecule name. Using the defaults from molecule1...")
                top12 += [(sec1, list1)]
            elif sec1 == "atoms":
                assert len(list1) == len(list2)
                list_temp = []
                for line1, line2 in zip(list1, list2):
                    list_temp += [line1 + [line2[1]] + line2[-2:]]
                top12 += [(sec1, list_temp)]
            elif sec1 == "bonds":
                bonds_1 = {(min(int(ls[0]), int(ls[1])), max(int(ls[0]), int(ls[1])), int(ls[2])): ls[3:]
                           for ls in list1}
                bonds_2 = {(min(int(ls[0]), int(ls[1])), max(int(ls[0]), int(ls[1])), int(ls[2])): ls[3:]
                           for ls in list2}
                bonds_1.update({key: val for key, val in bonds_2.items() if key not in bonds_1.keys()})
                bonds_2.update({key: val for key, val in bonds_1.items() if key not in bonds_2.keys()})
                top12 += [(sec1, [list(key) + val + bonds_2[key] for key, val in sorted(bonds_1.items())])]
            elif sec1 == "pairs":
                pairs = sorted(list({tuple(tuple(int(y) for y in x)) for x in list1} |
                                    {tuple(tuple(int(y) for y in x)) for x in list2}))
                top12 += [(sec1, pairs)]
            elif sec1 == "angles":
                angles_1 = {(min(int(ls[0]), int(ls[2])), int(ls[1]), max(int(ls[0]), int(ls[2])), int(ls[3])): ls[4:]
                            for ls in list1}
                angles_2 = {(min(int(ls[0]), int(ls[2])), int(ls[1]), max(int(ls[0]), int(ls[2])), int(ls[3])): ls[4:]
                            for ls in list2}
                angles_1.update({key: val for key, val in angles_2.items() if key not in angles_1.keys()})
                angles_2.update({key: val for key, val in angles_1.items() if key not in angles_2.keys()})
                top12 += [(sec1, [list(key) + val + angles_2[key] for key, val in sorted(angles_1.items())])]
            elif sec1 == "dihedrals":
                dihedr_1 = {tuple(int(x) for x in ls[:5]): ls[5:] for ls in list1}
                dihedr_2 = {tuple(int(x) for x in ls[:5]): ls[5:] for ls in list2}
                dihedr_1.update({key: [0, 0, 0] for key, val in dihedr_2.items() if key not in dihedr_1.keys()})
                dihedr_2.update({key: [0, 0, 0] for key, val in dihedr_1.items() if key not in dihedr_2.keys()})
                top12 += [(sec1, [list(key) + val + dihedr_2[key] for key, val in sorted(dihedr_1.items())])]
            elif sec1 == "system":
                print("Different system name. Using the defaults from system1...")
                top12 += [(sec1, list1)]
            elif sec1 == "molecules":
                print("Different molecule names / numbers. Using the defaults from system1...")
                top12 += [(sec1, list1)]

        with open(filename, "w") as file:
            for section, ls in top12:
                file.write("[ %s ]\n" % section)
                for line in ls:
                    file.write("\t" + "\t\t".join(map(str, line)) + "\n")
                file.write("\n")