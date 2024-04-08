import sys
import rdkit
from rdkit import Chem


def enumerateTorsions(mol):
   torsionSmarts = '[!$(*#*)&!D1]~[!$(*#*)&!D1]'
   torsionQuery = Chem.MolFromSmarts(torsionSmarts)
   matches = mol.GetSubstructMatches(torsionQuery)
   torsionList = []
   for match in matches:
     idx2 = match[0]
     idx3 = match[1]
     bond = mol.GetBondBetweenAtoms(idx2, idx3)
     jAtom = mol.GetAtomWithIdx(idx2)
     kAtom = mol.GetAtomWithIdx(idx3)
     if (((jAtom.GetHybridization() != Chem.HybridizationType.SP2)
       and (jAtom.GetHybridization() != Chem.HybridizationType.SP3))
       or ((kAtom.GetHybridization() != Chem.HybridizationType.SP2)
       and (kAtom.GetHybridization() != Chem.HybridizationType.SP3))):
       continue
     for b1 in jAtom.GetBonds():
       if (b1.GetIdx() == bond.GetIdx()):
         continue
       idx1 = b1.GetOtherAtomIdx(idx2)
       for b2 in kAtom.GetBonds():
         if ((b2.GetIdx() == bond.GetIdx())
           or (b2.GetIdx() == b1.GetIdx())):
           continue
         idx4 = b2.GetOtherAtomIdx(idx3)
         # skip 3-membered rings
         if (idx4 == idx1):
           continue
         torsionList.append((idx1, idx2, idx3, idx4))
   return torsionList

if (__name__ == "__main__"):
   mol = Chem.MolFromSmiles('CCCCCC')
   mol = Chem.AddHs(mol)
   print(mol)
   torsionList = enumerateTorsions(mol)
   for torsion in torsionList:
     sys.stdout.write(str(torsion) + '\n')
