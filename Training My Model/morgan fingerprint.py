import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

reactant_fps = []
mol = Chem.MolFromSmiles('[OH:101][C:1](=[O:2])[c:3]1[cH:4][c:5]([N+:6](=[O:7])[OH:8])[c:9]([S:10][c:11]2[c:12]([Cl:13])[cH:14][n:15][cH:16][c:17]2[Cl:18])[s:19]1')
if mol:
    print('success')