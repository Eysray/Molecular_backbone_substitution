from typing import List, Dict, Tuple, Optional # 导入类型提示
import pandas as pd
from rdkit import Chem # 导入环境

def atom_information(mol: Chem.Mol) -> List[Dict]:
    """获取分子中每个原子的信息"""
    # 获取分子中的原子
    atoms = mol.GetAtoms()
    atom_info = []
    for atom in atoms:
        info = {
        'symbol': atom.GetSymbol(),
        'atomicnum': atom.GetAtomicNum(),
        'index': atom.GetIdx(),
        'degree': atom.GetDegree(),
        'explicit_valence': atom.GetExplicitValence(),
        'implicit_valence': atom.GetImplicitValence(),
        'formal_charge': atom.GetFormalCharge(), #获取原子的形式电荷
        'total_num_hs': atom.GetTotalNumHs(),
        'total_degree': atom.GetTotalDegree(),
        'is_aromatic': atom.GetIsAromatic(),
        'isotope': atom.GetIsotope(),
        'mass': atom.GetMass(),
        'neighbors': atom.GetNeighbors(), #获取原子的邻居
        }
        atom_info.append(info)
    return atom_info

def bond_information(mol: Chem.Mol) -> List[Dict]:
    """获取分子中所有键的信息"""
    bonds = mol.GetBonds() # 获取所有键
    bond_info = []
    for bond in bonds:
        info = {
            'begin_atom_idx': bond.GetBeginAtomIdx(),
            'end_atom_idx': bond.GetEndAtomIdx(),
            'bond_type': bond.GetBondType().name,  # 使用枚举名称作为字符串
            'bond_idx': bond.GetIdx(),
        }
        bond_info.append(info)
    return bond_info

if __name__ == "__main__":
    smi = 'FC1=CC=C(C=C1[C@@H]2N(CCC2)C3=NC4=C(C=NN4C=C3)NC(OC5=CC=CC=C5)=O)F'
    mol = Chem.AddHs(Chem.MolFromSmiles(smi)) # 添加氢
    atom_info = atom_information(mol)
    df_atom = pd.DataFrame(atom_info)
    print(df_atom)