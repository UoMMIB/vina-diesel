# vina-diesel - an automated `autodock-vina` wrapper

[`Autodock vina`](http://vina.scripps.edu/) is a [molecular docking](https://en.wikipedia.org/wiki/Docking_(molecular)) command line program. **vina diesel** is a simple python wrapper that automates some of it.

## [original paper](https://pubmed.ncbi.nlm.nih.gov/19499576/)

```@article{trott2010autodock,
  title={AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading},
  author={Trott, Oleg and Olson, Arthur J},
  journal={Journal of computational chemistry},
  volume={31},
  number={2},
  pages={455--461},
  year={2010},
  publisher={Wiley Online Library}
}
```

## how `vina-diesel` works
`vina` works on `pdbqt` files, which `vina-diesel` creates in a tempfile using `openbabel`. `vina-diesel` requires a `pdb` file (receptor) and a smiles string (ligand).
Unless otherwise specified, all non-protein atoms are removed from the receptor structure when a `vdsl.protein` object is initialized (The first chain is also selected by default).

`vina` requires a box to be drawn around the area of the receptor to be docked. In `vina-diesel` this is specified using the residue numbers of the target site which can be found by playing with a 3D structure viewer like pymol. If the residue numbering in a structure is different to how you like (e.g. you're using a truncated or extended sequence) then provide the full sequence to `vdsl.protein` to align that sequence to the structure and you'll be able to use your sequence numbering.

`vina-diesel` calls the `vina` binary as a subprocess and stores the results in a `results` object, which can be used to save the poses and energies or access a `pandas` DataFrame of the pose energies and a `biopandas.pdb` object containing the coordinates. The results object is a little bit scrappy, sorry about that.

## simple example

```python
import vdsl

p = vdsl.protein('3b4y.pdb') ###
results = p.dock('[O-][N+](=O)c1cn2C[C@@H](COc2n1)OCc3ccc(OC(F)(F)F)cc3', 
    target_sites = [38, 44, 175, 176, 199, 252, 256, 261, 283, 311]) # 

results.save('thatwaseasy')
```

## another example


```python
import vdsl

MY_SEQUENCE = 'MTAJSKNFHASBOYBRWBRHA...' # keyboard mash
COFACTORS = ['HEM', 'FAD'] # cofactors to keep -- names as they appear in the pdb file

p = vdsl.protein('3b4y.pdb',
		seq = MY_SEQUENCE,
		keep = COFACTORS)
results = p.dock('[O-][N+](=O)c1cn2C[C@@H](COc2n1)OCc3ccc(OC(F)(F)F)cc3', 
	    target_sites = [38, 44, 175, 176, 199, 252, 256, 261, 283, 311], # 
		exhaustiveness=16)  # integer 1-16 (coarse-fine grained)

results.save('thatwaseasy')

pritn(results.scores) # pd.DataFrame
```
