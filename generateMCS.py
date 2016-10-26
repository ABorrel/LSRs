from nams import nams
from rdkit import Chem
import openbabel
import os

from fmcs import fmcs
import loader
import runSoftwares


def BuildMatrixSimilarity(d_complex_ref, d_complex_pose, type_MCS, p_filout):
    # test if file exist

    if os.path.exists(p_filout):
        print
        "INNNN", os.path.getsize(p_filout)
        if os.path.getsize(p_filout) > 50:
            return loader.LoadMCSResult(p_filout)

    d_out = {}
    for lig_pose in d_complex_pose.keys():
        pose_pose = d_complex_pose[lig_pose].keys()[0]
        for lig_ref in d_complex_ref.keys():
            pose_ref = d_complex_ref[lig_ref].keys()[0]

            print
            lig_ref, pose_ref, lig_pose, pose_pose

            if type_MCS == "NAMS":
                similarity = SimilarityNAMSJaccard(d_complex_ref[lig_ref][pose_ref]["SMILE CODE"],
                                                   d_complex_pose[lig_pose][pose_pose]["SMILE CODE"])
            if type_MCS == "FMCS":
                similarity = SimilarityFMCSTanimoto(d_complex_ref[lig_ref][pose_ref]["SMILE CODE"],
                                                    d_complex_pose[lig_pose][pose_pose]["SMILE CODE"])

            if not lig_ref in d_out.keys():
                d_out[lig_ref] = {}
            d_out[lig_ref][lig_pose] = similarity

    p_matrix = WriteMatrixSimilarity(d_out, p_filout)

    # run color matrix
    runSoftwares.CardMatrix(p_matrix)

    # append matrix image => todo

    return d_out


def SimilarityNAMSJaccard(smile1, smile2):
    ms = nams.Nams()
    mol_t1 = ("smi", smile1)
    mol_t2 = ("smi", smile2)
    mol1, mol_info1 = ms.get_mol_info(mol_t1[0], mol_t1[1])
    mol2, mol_info2 = ms.get_mol_info(mol_t2[0], mol_t2[1])

    # similarity combination
    sim12, d_atoms = ms.get_similarity(mol_info1, mol_info2)
    sim11, d_atoms = ms.get_similarity(mol_info1, mol_info1)
    sim22, d_atoms = ms.get_similarity(mol_info2, mol_info2)

    # based on a Jaccard score
    score_Jaccard = sim12 / (sim11 + sim22 - sim12)

    return score_Jaccard


def SimilarityFMCSTanimoto(smile1, smile2):
    mols = [Chem.MolFromSmiles(smile1), Chem.MolFromSmiles(smile2)]
    sim12 = fmcs(mols).num_atoms

    nb1 = mols[0].GetNumAtoms()
    nb2 = mols[1].GetNumAtoms()

    # based on tanimoto score on atom commun
    score_tanimoto = float(sim12) / float((nb2 - sim12) + (nb1 - sim12) + sim12)

    return score_tanimoto


def WriteMatrixSimilarity(d_in, p_filout):
    filout = open(p_filout, "w")

    l_lig_ref = d_in.keys()
    l_lig_test = d_in[d_in.keys()[0]].keys()
    filout.write("\t".join(d_in.keys()) + "\n")

    for lig_test in l_lig_test:
        filout.write(str(lig_test))
        for lig_ref in l_lig_ref:
            filout.write("\t" + str(d_in[lig_ref][lig_test]))
        filout.write("\n")
    filout.close()

    return p_filout