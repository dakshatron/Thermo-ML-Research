# Problem Statement
- the NIST Reaction Kinetics database has a treasure trove of data on how quickly many different reactions progress
- a new extremely high quality dataset on reaction pathways of organic reactions has 30K+ reaction mechanisms/steps broken down by hand by chemists, in standardized (canonical SMILES) formats, specifying exactly which bonds are broken/formed, assigning IDs for you to follow a single atom through a reaction

I'm trying to leverage both to predict a given reaction's Arrhenius rate constant equation parameters (k = A * exp(- E_a/RT, I'm predicting A and E_a, the pre-exponential factor and activation energy, respectively) by:
- decomposing any given reaction into the likely mechanisms/steps in between
- homing in on the specific atoms most likely to speed up, or slow down, a reaction
- getting an idea of how fast or slow certain mechanisms are.

The database and publication links are here:
- https://kinetics.nist.gov/kinetics/
- https://www.nature.com/articles/s41597-024-03709-y


# Current Progress
I'm still playing around with the 'Pathways Dataset', currently with scikit (Random Forest to be specific) to attempt to guess which molecules are the center of the reaction using Morgan fingerprints (2048 length bit (0/1) vector 'simulating' a circular radius around a 2D map (similar to nodes of a graph neural network) by hashing and 'folding)

I'm also experimenting with the MACCS datatype, which's a 166 bit vector, as Morgan fingerprints lead to very high dimensional, and as a result very sparse vector spaces

I'm having trouble integrating my code with the MechFinder github repo I cloned from the Nature publication's backend, will try and fix/push it by the weekend, or just make a new repo, though I want to keep it all in one VSCode workspace


