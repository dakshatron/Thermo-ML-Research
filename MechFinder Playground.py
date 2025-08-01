from MechFinder import MechFinder
import numpy as np
import pandas
from rdkit import Chem


finder = MechFinder(collection_dir='collections')
sampled_rxns = pandas.read_csv('data/samples.csv')['reaction'].tolist() # randomly sampled reactions
sampled_rxn = np.random.choice(sampled_rxns, 1)[0]
Chem.Draw.MolsToGridImage([Chem.MolFromSmiles(smi) for smi in  sampled_rxn.split('>>')], subImgSize=(400, 400))

