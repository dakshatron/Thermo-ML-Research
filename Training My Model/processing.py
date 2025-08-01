import pandas
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

from typing import List, Dict, Any, Optional, Tuple # for specifying data types so things dont get fucked up

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
from rdkit.Chem import rdChemReactions
from rdkit.DataStructs import ConvertToNumpyArray

pathwaysDataframe = pandas.read_csv('Pathways Dataset.csv')
processedDataframe = pandas.DataFrame()

# splitting reactions into seperate reactants and products
fingerprintsArray: list[np.ndarray] = []
processedRows: List = []
for rowIndex, rowContents in pathwaysDataframe.iterrows():
reaction = rowContents['Reaction Mechanism']
# automatically parses (splits rxn based on >> and .) and converts to "rxn object"
rxnObj = rdChemReactions.ReactionFromSmarts(reaction, useSmiles=True)

# defines reaction fingerprint parameters object as that's all CreateDifferenceFingerprintForReaction takes
rxnFingerprintParams = rdChemReactions.ReactionFingerprintParams()
rxnFingerprintParams.fpSize = 4096
rxnFingerprintParams.fpType = AllChem.FingerprintType.MorganFP # captital FP, doesn't recognize otherwise!!!
# rxnFingerprintParams.fpRadius = 2

morganDiffVector = rdChemReactions.CreateDifferenceFingerprintForReaction(rxnObj, rxnFingerprintParams) # takes rxnObj, creates difference vector
morganDiffNpArray = np.zeros((4096,), dtype=np.int32) # empty placeholder array
ConvertToNumpyArray(morganDiffVector, morganDiffNpArray) # fills morganDiffNpArray with morganDiffVector data

fingerprintsArray.append(morganDiffNpArray) # appends the fingerprint to the list
processedRows.append(rowContents) # appends the row contents to the list

# --- Assemble the Final DataFrame (The Efficient Way) ---
print(f"\nSuccessfully processed {len(processedRows)} reactions.")

# 1. Create the fingerprint DataFrame from your list of arrays
fingerprintDataframe = pandas.DataFrame(
fingerprintsArray,
columns=[f'bit_{i}' for i in range(4096)],
index=processedRows # Use the saved indices to align data correctly
)

# 2. Get the original data that corresponds to the successful fingerprints
original_data_succeeded = pathwaysDataframe.loc[processedRows]

# 3. Concatenate the original data with the new fingerprint features
finalDataframe = pandas.concat([original_data_succeeded, fingerprintDataframe], axis=1)

finalDataframe.to_csv('Morgan Difference Pathways Dataset.csv', index=False)