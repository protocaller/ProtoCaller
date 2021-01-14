import logging
import os
logging.basicConfig(level=logging.INFO)

from ProtoCaller.Utils.fileio import Dir
from ProtoCaller.Ensemble import Ligand, Protein, Ensemble, Perturbation

from glob import glob
from BioSimSpace.IO  import saveMolecules
from ProtoCaller.Wrappers.babelwrapper import babelTransform
import ProtoCaller.Utils.runexternal as _runexternal

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdMolAlign import AlignMol

def generate_teth(ref_mol, target_mol):  # generate sd file for tethred docking
    reference = Chem.MolFromMolFile(ref_mol, removeHs=True)
    ligands = Chem.SDMolSupplier(target_mol,removeHs=True)
    w = Chem.SDWriter('teth.sd')
    for mol in ligands:
        mols=[reference,mol]
        mcsResult=rdFMCS.FindMCS(mols,threshold=0.9,completeRingsOnly=True) #find the maximum common substructure
        
        if mcsResult.smartsString and len(mcsResult.smartsString)>0 :
            patt = Chem.MolFromSmarts(mcsResult.smartsString,mergeHs=True)
            ref = AllChem.ReplaceSidechains(reference, patt)
            if ref:
                core=AllChem.DeleteSubstructs(ref,Chem.MolFromSmiles('*'))
                core.UpdatePropertyCache()
                newmol=Chem.Mol(mol)	#create a new instance of the molecule, as we will change the coordinates
                AllChem.ConstrainedEmbed(newmol,core)	#constrained minimization of newmol versus the core of the reference
                rms=float(newmol.GetProp('EmbedRMS'))
                print('Common core: ',Chem.MolToSmiles(core),'  - RMSD:',rms)
                tethered_atom_ids=newmol.GetSubstructMatches(patt)	#that's to get the atom ids only
                if tethered_atom_ids:
                    t=tethered_atom_ids[0]
                    t1=map(lambda x:x+1, list(t))
                    ta=','.join(str(el) for el in t1)
                    nm=Chem.AddHs(newmol, addCoords=True)	#create a new 3D molecule  and add the TETHERED ATOMS property
                    nm.SetProp('TETHERED ATOMS',ta)
                    w.write(nm)
                    w.close()

def generate_prm(receptor_file, ref_mol):  # generate parameter file for rdock
    with open('../template.prm', 'rt') as f:
        content = f.read()
    new_content = content.format(receptor_file=receptor_file, ref_mol=ref_mol)
    with open('cavity.prm', 'wt') as f:
        f.write(new_content)

def prepare_system(ligands) -> None: 
    with Dir('SYSTEM', overwrite=True):
        pdb_id = '3eyg'
        protein_proto = Protein(pdb_id, ligand_ref='1', workdir='before_docking')

        ligands_proto = []
        for idx, item in enumerate(ligands):
            li = Ligand(item, protonated=False, minimise=False, name=f"ligand{idx}", workdir="Ligands")
            ligands_proto.append(li)
            ligands_proto[-1].protonate()
            ligands_proto[-1].parametrise()     # here we got the parametrised sdf file for docking

        protein_proto.filter(ligands=None, waters="all")
        protein_proto.prepare()
        protein_proto.parametrise()

        ligand_ref = glob(f'before_docking/*.sdf')[0]
        saveMolecules('complex', protein_proto.complex_template, 'pdb')   # save the complex as pdb file
        babelTransform('complex.pdb')  # converted to mol2 file required by rdock
        generate_prm('complex.mol2', ligand_ref)  # generate parameter file needed by rdock
        _runexternal.runExternal('rbcavity -r cavity.prm -W -d > generate_cavity.log')  # generate cavity
        generate_teth(ligand_ref, 'Ligands/ligand0_protonated.sdf')  # generate tethred docking configuration in sd file
        _runexternal.runExternal(f'rbdock -r cavity.prm -p dock.prm -n 100 -i teth.sd -o docking_out > docking_out.log ')  # docking
        _runexternal.runExternal("sdsort -n -s -fSCORE docking_out.sd | sdfilter -f'$_COUNT == 1' > result.sd") # sort the result and get the best pose

        os.rename('result.sd', 'result.sdf')   # .sdf was needed by protocaller
        docked_ref = Ligand('result.sdf')
        docked_protein = Protein(pdb_id, ligand_ref=docked_ref, workdir='after_docking')   # reinitialize an protein object

        perts = []
        pert = Perturbation(ligands_proto[0], ligands_proto[1])
        pert.alignToReference(docked_protein.ligand_ref)#, mcs_parameters=dict(break_recursively=False))
        pert.alignToEachOther()#mcs_parameters=dict(keep_stereo=False))
        perts += [pert]

        docked_system = Ensemble("GROMACS", protein=docked_protein, morphs=perts,
                          box_length_complex=7, ligand_ff="gaff2",
                          workdir=docked_protein.workdir.path)

        docked_system.protein.filter(ligands=None, waters="all")
        docked_system.protein.prepare()
        docked_system.protein.parametrise()

        docked_system.prepareComplexes(intermediate_files=True)   # save intermediate file for diagnose

if __name__ == '__main__':
    ligands = ['CC(C)[C@H](CNC(=O)CC#N)NC1=C2C=CNC2=NC=N1',
            'O=C(CC#N)NCCNC1=C2C=CNC2=NC=N1']

    prepare_system(ligands)
