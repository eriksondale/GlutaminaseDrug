# Written by Dale Erikson

# The purpose of this script is to determine what molecules a library are
# able to be formed via a given rxn

# note: Rxns are expected to be given in the reverse direction
import glob
from sys import argv as argument
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import rdBase



def validRxn(reactant, reaction):
	try:
		product = None
		product = reaction.RunReactants(reactant)
		if(len(product) == 0):
			#print('Product length is 0: ' + Chem.MolToSmiles(reactant[0]) + ' ' + reactionPlan)
			return False
		else:
	#	return False

        #Check that forward and reversal are same
        #compound = Chem.MolToSmiles(reactant[0])
			reverseReactionPlan = reactionPlan[reactionPlan.find('>>')+2:] + '>>' + reactionPlan[0:reactionPlan.find('>>')]
			reverseReaction = AllChem.ReactionFromSmarts(reverseReactionPlan)

			for prods in product:
				try:
					for prod in prods:
						Chem.SanitizeMol(prod)
					reactant = reverseReaction.RunReactants(prods)
					for pairs in reactant:
						for molecule in pairs:
							moleculeBit = FingerprintMols.FingerprintMol(molecule)
							compoundBit= FingerprintMols.FingerprintMol(mol)
							similarity = DataStructs.FingerprintSimilarity(moleculeBit, compoundBit)
							#print(similarity)
							if(similarity == 1): # Rxn is valid b/c product and reverse product is found to be same
								return True
				except:
					pass

		   	#print("Failure to return same molecule: " + Chem.MolToSmiles(reactant[0]) + " " + reactionPlan)
		  	return False
	except:
		#print('Overall error')
		return False


if(len(argument) != 3):
	print("Error in argument length; molecule and rxn files needed....")
else:
	if(argument[1].endswith('.sdf')):
		suppl = Chem.SDMolSupplier(argument[1])
		molList = [x for x in suppl if x is not None]
	elif(argument[1].endswith('.smi')):
		molList = []
		tempFile = open(argument[1],"r")
		for line in tempFile:
			line = line.strip("\n")
			if line is not None:
				molList.append(Chem.MolFromSmiles(line))
	elif(argument[1].endswith('.pdb')):
		molList = []
		molList.append(Chem.MolFromPDBFile(argument[1]))
	else:
		print('Need .sdf or .smi or .pdb file to read')
		sys.exit()
	with open(argument[2],"r") as RxnFile:
		for line in RxnFile:
			numInvalid = 0
			numValid = 0
			line = line.split("\t")
			reactionPlan = line[1]
			reaction = AllChem.ReactionFromSmarts(reactionPlan)
			for mol in molList:
				if mol is None:
					pass
				else:
					Chem.SanitizeMol(mol)
					reactant = []
					reactant.append(mol)
					validity = validRxn(reactant, reaction)
					if validity:
						numValid += 1
						print(Chem.MolToSmiles(mol))
					else:
						numInvalid += 1
			print(line[0] + " valid rxns: " + str(numValid))
			print(line[0] + " invalid rxns: " + str(numInvalid))
