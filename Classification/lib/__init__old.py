from key_generation import KeyGeneration
from pairwiseLocalAlignment import FeatureSelection
from vectFromLocalFeatures import Vectorization
#Uncomment the line below for sequence execution - not parallel
from jacOnlyFromVect import JaccardCoefficient
#Uncomment the line below for parallel
#from jacOnlyFromVect_parallel import JaccardCoefficient
#no changes below
from dendo_old import Dendograming
from similarityScore import SimilarityScore