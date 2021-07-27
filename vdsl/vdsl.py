import tempfile
import os
import shutil
from distutils.spawn import find_executable
import subprocess
import re
from itertools import chain

import pandas as pd

from biopandas.pdb import PandasPdb
import nwalign3 as nw
from openbabel import pybel

pybel.ob.obErrorLog.SetOutputLevel(0)

class mol:
    '''
    protein & results poses inherit from this class
    '''
    def __init__(self, struc):
        self.struc = struc
    @property
    def df(self):
        data = PandasPdb().read_pdb(self.struc)
        return data.df['ATOM'].append(data.df['HETATM']) # all atoms
    def save(self, save_path):
        shutil.copyfile(self.struc, save_path)


class protein(mol):
    def __init__(self, struc, seq = None, keep = [], key_sites = []):
        super().__init__(struc)
        self.key_sites = key_sites
        self.keep = keep
        self.CACHE = tempfile.mkdtemp()
        self.struc = pdb_fns.clean_pdb(struc = self.struc,
                    save_path = os.path.join(self.CACHE, 'clean.pdb'),
                    keep = self.keep)
        self.pdb_seq = pdb_fns.get_seq(self.struc)
        self.seq = self.pdb_seq if seq == None else seq
    
    def dock(self, 
            smiles, 
            save_path = None,
            target_sites = None,
            exhaustiveness = 8):
        if target_sites is None:
            target_sites = self.key_sites
        results = vina.dock(self.struc,
                    smiles,
                    save_path = save_path,
                    keep = self.keep,
                    target_sites = target_sites,
                    exhaustiveness = exhaustiveness)
        return results

class pdb_fns:
    def clean_pdb(struc, save_path, keep = [], chain_selection = None):
        structure = PandasPdb().read_pdb(struc)
        atoms = structure.df['ATOM'].copy()

        # added 20210423 - there were chain selection problems
        chains = atoms['chain_id'].unique()
        if chain_selection is None:
            chain_selection = chains[0]
        #### 

        hetatms = structure.df['HETATM'].copy()
        atoms = atoms.loc[atoms['chain_id'] == chain_selection,:]
        # edge case - het atoms didnt have chain assigned (single chain struc - 3iw2)
        if chain_selection in hetatms['chain_id']:
            hetatms = hetatms.loc[hetatms['chain_id'] == chain_selection,:]
        het_garbage = [i for i in hetatms['residue_name'].unique() if i not in keep]
        hetatms = hetatms.loc[hetatms['residue_name'].isin(het_garbage) == False,:]
        structure.df['ATOM'] = atoms
        structure.df['HETATM'] = hetatms
        structure.to_pdb(save_path)
        return save_path

    def get_seq(struc):
        structure = PandasPdb().read_pdb(struc)
        sequences = structure.amino3to1() # cols = ['chain_id', 'residue_name']
        seqs = [''.join(sequences.loc[sequences['chain_id'] == i,'residue_name'].to_list()) for i in sequences['chain_id'].unique()]
        
        return seqs[0] if len(seqs) == 1 else seqs

    def draw_box(struc, key_sites):
        receptor = PandasPdb().read_pdb(struc)
        df = receptor.df['ATOM']
        target_site = df.loc[df['residue_number'].isin(key_sites),:]
        coords = target_site.loc[:,['x_coord','y_coord','z_coord']]
        center = coords.mean(axis=0)
        sizes = (coords.max(axis=0) - coords.min(axis=0)) * 1.2
        box = {'--center_x':center['x_coord'],
                '--center_y':center['y_coord'],
                '--center_z':center['z_coord'],
                '--size_x':sizes['x_coord'],
                '--size_y':sizes['y_coord'],
                '--size_z':sizes['z_coord']}
        return box

class obabel_fns:
    def pdb_to_pdbqt(pdb, save_path):
        m = list(pybel.readfile('pdb',pdb))
        assert len(m) == 1
        m = m[0]
        m.addh()
        m.write('pdbqt', save_path, opt={'r':True}, overwrite=True) # opt:r = rigid - less errors?? - revisit this
        return save_path

    def smiles_to_pdbqt(smiles, save_path):
        m = pybel.readstring('smi',smiles)
        m.OBMol.StripSalts()
        m.addh()
        m.make3D()
        m.write('pdbqt',save_path, overwrite=True)
        return save_path

    def pdbqt_to_pdb(pdbqt, save_path):
        m = list(pybel.readfile('pdbqt',pdbqt))
        assert len(m) == 1
        m = m[0]
        # already 3d
        m.write('pdb', save_path, overwrite=True)


class vina:
    def dock(receptor_pdb,
            smiles,
            save_path = None,
            keep = [],
            target_sites = [],
            exhaustiveness=8,
            vina_executable = find_executable('vina'),
            vina_split_executable = find_executable('vina_split')):
        # check there's a box
        if target_sites == []:
            raise Exception('no target residues selected')

        CACHE = tempfile.mkdtemp()
        raw_vina_results = os.path.join(CACHE, 'vina.result')
        # todo : if not clean  - check if dock from protein object
        # otherwise this executes twice
        clean_receptor_pdb = pdb_fns.clean_pdb(receptor_pdb, os.path.join(CACHE, f'{os.path.basename(receptor_pdb)}.clean'), keep = keep)
        receptor_pdbqt = obabel_fns.pdb_to_pdbqt(clean_receptor_pdb, os.path.join(CACHE,'receptor.pdbqt'))
        ligand_pdbqt = obabel_fns.smiles_to_pdbqt(smiles, os.path.join(CACHE,'ligand.pdbqt'))
        args = {'--receptor':receptor_pdbqt,
                    '--ligand':ligand_pdbqt,
                    '--out':raw_vina_results,
                    '--exhaustiveness':exhaustiveness}
        
        args.update(pdb_fns.draw_box(clean_receptor_pdb, target_sites)) # add box dims to args

        args_list_vina = [vina_executable] + [str(i) for i in chain.from_iterable(args.items())]

        # execute
        p1 = subprocess.check_output(args_list_vina)
        
        # create results object
        # clean_receptor_pdb to pdb
        docking_scores = vina.extract_scores(p1.decode())
        poses = vina.vina_split(raw_vina_results, vina_split_executable)
        results = vina.results(clean_receptor_pdb, [os.path.join(poses, i) for i in os.listdir(poses)], docking_scores)
        
        return results

    def vina_split(raw_vina_results, vina_split_executable):
        # vina_split
        args_list_vina_split = [vina_split_executable, '--input', raw_vina_results]
        p = subprocess.Popen(args_list_vina_split, stdout=subprocess.DEVNULL)
        p.wait() # outputs odbqt files
        results_dir = os.path.dirname(raw_vina_results)
        poses = [os.path.join(results_dir, i) for i in os.listdir(results_dir) if 'vina.result_ligand' in i]
        clean_results = os.path.join(results_dir, 'pose_pdbs')
        os.mkdir(clean_results)
        for i in poses:
            save_path = os.path.basename(i).replace('pdbqt','pdb')
            obabel_fns.pdbqt_to_pdb(i, os.path.join(clean_results, save_path))
        return clean_results # path to foder containing pdb of poses

    def extract_scores(text):
        # extract scores from vina output
        text = text.split('\n')
        table_start = ['---+--' in i for i in text].index(True) + 1
        is_all_ints = lambda l : sum([re.search('-?\d+', i) is not None for i in l])
        table = []
        for row in text[table_start:]:
            items = row.split()
            if len(items) == 4 and is_all_ints(items):
                table.append(dict(zip(['mode','affinity (kcal/mol)', 'dist from best mode - rmsd - ub','dist from best mode - lb'], items)))
        return pd.DataFrame(table)

    class results:
        '''
        poses & score df
        '''
        def __init__(self, receptor, poses, scores):
            self.poses = {int(re.findall('\d+',os.path.basename(i))[0])\
                    :mol(i) for i in poses}
            self.receptor = receptor # path to clean pdb
            self.scores = scores.astype(float)
            self.dictionary = {os.path.basename(i.struc):{'mol':i, 'affinity':j} for i,j in zip(self.poses.values(), self.scores['affinity (kcal/mol)'])}
        def save(self, save_path):
            os.makedirs(save_path, exist_ok = True)
            self.scores.to_csv(os.path.join(save_path, 'scores.csv'))
            for i in self.poses:
                pose_i = self.poses[i]
                pose_i.save(os.path.join(save_path, os.path.basename(pose_i.struc)))
            # saves pdb
            shutil.copyfile(self.receptor, os.path.join(save_path, 'clean_receptor.pdb'))



class utils:
    def aln(s1, s2):
        aln1, aln2 = nw.global_align(s1,s2)
        return aln1, aln2

    def diff(s1,s2):
        # s1 - canonical seq
        # s2 - pdb seq
        # offset needed to map aligned positions to pyrosetta positions
        # todo: test where len(s1) < len(s2)
        offset = lambda s, idx : sum([i == '-' for i in s[:idx]])
        return {i - offset(s2,i):{'from':x, 'to':y} for i, (x,y) in enumerate(zip(s2,s1) ,1) if x != y and x != '-' and y != '-'}


