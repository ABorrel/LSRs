from nams import nams



def get_Tanimoto(psmile1, psmile2):
    ms = nams.Nams()

    fsmile1 = open(psmile1, "r")
    smile1 = fsmile1.readlines()[0].strip()
    fsmile1.close()

    fsmile2 = open(psmile2, "r")
    smile2 = fsmile2.readlines()[0].strip()
    fsmile2.close()


    print(smile1, smile2)
    mol_t1 = ("smi", smile1)
    mol_t2 = ("smi", smile2)
    mol1, mol_info1 = ms.get_mol_info(mol_t1[0], mol_t1[1])
    mol2, mol_info2 = ms.get_mol_info(mol_t2[0], mol_t2[1])

    # similarity combination
    sim12, d_atoms = ms.get_similarity(mol_info1, mol_info2)
    sim21, d_atoms = ms.get_similarity(mol_info2, mol_info1)
    sim11, d_atoms = ms.get_similarity(mol_info1, mol_info1)
    sim22, d_atoms = ms.get_similarity(mol_info2, mol_info2)

    # based on a Jaccard score
    score12 = sim12 / (sim11 + sim22 - sim12)
    score21 = sim21 / (sim11 + sim22 - sim21)

    return [score12, score21]
