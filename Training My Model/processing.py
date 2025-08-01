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
columns_list = ['Reactant 1', 'Reactant 2', 'Product 1', 'Product 2', 'Mechanism Class', 'Mechanism Bond Changes', 'Reaction Mechanism']
processedDataframe = pandas.DataFrame(columns=columns_list)

# splitting reactions into seperate reactants and products

def smiles2morgan(smilesParam: str, radiusParam: Optional[int] = 2, nBitsParam: Optional[int] = 2048) -> np.ndarray:
    """
    converts a SMILES string to a Morgan fingerprint bit (ie. 1s and 0s) vector 2048 entries long

    Args:
        smilesParam (str): the SMILES string to convert
        radiusParam (int, optional): the radius of the Morgan fingerprint. Defaults to 2
        nBitsParam (int, optional): the length of the fingerprint vector. Defaults to 2048
    
    Returns:
        np.ndarray: a numpy array of 1s and 0s representing the Morgan fingerprint. Will save to csv something like this:
        [Row 1: [1, 0, 0, 1, ..., 'Mechanism Class', 'Mechanism Bond Changes']]
    """
    mol = Chem.MolFromSmiles(smilesParam)
    if mol is None:
        return np.array([])  # empty array if SMILES invalid
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radiusParam, nBits=nBitsParam)
    return np.array(fingerprint)

smarts2morganDiff = rdChemReactions.ReactionFromSmarts()

for rowIndex, rowContents in pathwaysDataframe.iterrows():
    reaction = rowContents['Reaction Mechanism']
    # automatically parses (splits rxn based on >> and .) and converts to "rxn object"
    smarts2morganDiff = rdChemReactions.ReactionFromSmarts()

    morganDiffVector = rdChemReactions.CreateDifferenceFingerprintForReaction(reaction, fpSize=4096, fpType=AllChem.MorganFp, radius=2)
    morganDiffNpArray = np.zeros((4096,), dtype=np.int32) # empty placeholder array
    ConvertToNumpyArray(morganDiffVector, morganDiffNpArray) # fills morganDiffNpArray with morganDiffVector data


    ##### Found RDKit function that does this for me (CreateDifferenceFingerprintForReaction) #####
    # rxnParts: list = reaction.split('>>')

    # if len(rxnParts) > 1:
    #     reactants:str = rxnParts[0]
    #     products:str = rxnParts[1]
        
    #     reactantList: list = [s.strip() for s in reactants.split('.')]
    #     productList: list = [s.strip() for s in products.split('.')]
        
    #     rowReactantFingerprints: List[np.ndarray] = []
    #     rowProductFingerprints: List[np.ndarray] = []

    #     for reactant in reactantList:
    #         reactantFingerprint = smiles2morgan(reactant, radiusParam=2, nBitsParam=2048)
    #         rowReactantFingerprints.append(reactantFingerprint)
        
    #     for product in productList:
    #         productFingerprint = smiles2morgan(product, radiusParam=2, nBitsParam=2048)
    #         rowReactantFingerprints.append(productFingerprint)


    #     processedDataframe.loc[rowIndex, 'Reactant 1'] = rowReactantFingerprints[0] if len(rowReactantFingerprints) > 0 else None
    #     processedDataframe.loc[rowIndex, 'Reactant 2'] = reactantList[1] if len(reactantList) > 1 else None
    #     processedDataframe.loc[rowIndex, 'Product 1'] = productList[0] if len(productList) > 0 else None
    #     processedDataframe.loc[rowIndex, 'Product 2'] = productList[1] if len(productList) > 1 else None

    #     processedDataframe.loc[rowIndex, 'Mechanism Class'] = rowContents['Mechanism Class']
    #     processedDataframe.loc[rowIndex, 'Mechanism Bond Changes'] = rowContents['Mechanism Bond Changes']
    #     processedDataframe.loc[rowIndex, 'Reaction Mechanism'] = rowContents['Reaction Mechanism']
    # else:
    #     continue

processedDataframe.to_csv('Processed Pathways Dataset.csv', index=False)