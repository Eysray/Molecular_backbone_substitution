import re
from rdkit import Chem # 导入环境
import pandas as pd # 导入环境
import random
from typing import List, Dict, Tuple, Optional # 导入类型提示
from itertools import chain

'''
获取分子侧链数量
'''
def get_cut_and_match_information(mol: Chem.Mol, core: Chem.Mol) -> Optional[Dict[str, int]]:
    Chem.Kekulize(mol, clearAromaticFlags=True)
    if mol.HasSubstructMatch(core) == True:
        m_cut = Chem.ReplaceCore(mol, core)
        side_chains = Chem.GetMolFrags(m_cut, asMols = True, sanitizeFrags = True)
        output_dict = {'side_chains_num': len(side_chains), #侧脸数量，不对应接口数
                'SubstructMatch_num': len(mol.GetSubstructMatches(core)), #匹配数量，对应官能团数量
                }
        matchatoms_index = mol.GetSubstructMatches(core) #匹配原子索引
        flattened_list = list(chain.from_iterable(item if isinstance(item, tuple) else (item,) for item in matchatoms_index)) #展平匹配原子元组索引
        flattened_tuple = tuple(flattened_list) #所有匹配原子索引的元组
        return output_dict, flattened_tuple
    else:
        return None
    
'''
一次只替换一个官能团。
均接受SMILES码输入
mol_smi: 分子
core_sma: 母核
new_core_smi: 新母核
返回值:
result_smi: 替换后的分子SMILES码
result_mol: 替换后的分子Chem.Mol对象
注意：
1. new_core_smi需要在接口添加虚拟原子
2. 分子在切掉母核后产生的碎片数目要与新母核的接口数相同
'''
def replace_core(mol_smi: str, core_sma: str, new_core_smi: str, useChemMol = False) -> Tuple[str, Chem.Mol]:
    ## 准备分子
    if useChemMol: # 如果useChemMol为真，输入Chem.Mol对象
        mol = mol_smi
        core = core_sma
        new_core = new_core_smi
    else:
        mol = Chem.MolFromSmiles(mol_smi)
        core = Chem.MolFromSmarts(core_sma) # 母核需要是SMARTS格式，可以匹配更多结构
        new_core = Chem.MolFromSmiles(new_core_smi) #new_core_smi需要在接口添加虚拟原子，Isotope分别设置依次增加
    Chem.Kekulize(mol, clearAromaticFlags=True) # 清除芳香性,以匹配SMART匹配检索，注意：共振式无法检索。
    ## 切割分子，并把片段和新母核合并在一个对象里
    m_cut = Chem.ReplaceCore(mol, core) # 一次只替换一个官能团，一次只切割一个。
    side_chains = Chem.GetMolFrags(m_cut, asMols = True, sanitizeFrags = True)
    combined = new_core
    for chain in side_chains:
        combined = Chem.CombineMols(chain, combined)
    ## 根据原子属性标记出接口虚拟原子
    isotope_dict = {1: []} # {虚拟原子Isotope:虚拟原子Idx},使用同位素虚拟原子来指定接口
    for atom in combined.GetAtoms():
        isotope = atom.GetIsotope()
        symbol = atom.GetSymbol()
        if isotope in isotope_dict and symbol == '*':
            isotope_dict[isotope].append(atom.GetIdx())
        if isotope not in isotope_dict and symbol == '*':
            isotope_dict[isotope] = [atom.GetIdx()]
    ## 相同标记的虚拟原子之间成键
    rw_mol = Chem.RWMol(combined)
    for key in isotope_dict.keys():
        isotope_key = isotope_dict[key]
        rw_mol.AddBond(isotope_key[0], isotope_key[1], Chem.BondType.SINGLE)
    modified_mol = rw_mol.GetMol() # 获取修改后的分子
    Chem.SanitizeMol(modified_mol) # 检查并清理分子
    ## 正则匹配式删除虚拟原子
    modified_smi = Chem.MolToSmiles(modified_mol)
    pattern = r'\[(\d+)\*\]|\(\[(\d+)\*\]\)'
    result_smi = re.sub(pattern, '', modified_smi)
    result_mol = Chem.MolFromSmiles(result_smi)
    return result_smi, result_mol

def replace_core_randomvision(mol: Chem.Mol, datacsv_path: str) -> Tuple[str, Chem.Mol]:
    Chem.Kekulize(mol, clearAromaticFlags=True) # 清除芳香性,以匹配SMART匹配检索，注意：共振式无法检索。
    df = pd.read_csv(datacsv_path) 
    quarys = df['no_star'].values.tolist()
    match_list = []
    # 子结构检查
    for index, quary in enumerate(quarys):
        check = mol.HasSubstructMatch(Chem.MolFromSmarts(quary))
        if check == True:
            match_list.append(index) # 生成包含子结构的quarys的索引列表
    choice_index = random.choice(match_list)
    core = Chem.MolFromSmarts(quarys[choice_index])   
    ## 切割分子，并把片段和新母核合并在一个对象里
    m_cut = Chem.ReplaceCore(mol, core) # 一次只替换一个官能团，一次只切割一个。
    side_chains = Chem.GetMolFrags(m_cut, asMols = True, sanitizeFrags = True)
    side_chains_combnind = Chem.RWMol() #empty_mol 
    for chain in side_chains:
        side_chains_combnind = Chem.CombineMols(chain, side_chains_combnind)
    num_stars = Chem.MolToSmiles(side_chains_combnind).count('*') # 计算接口数量
    #根据接口数量挑选新母核
    new_coredf = df[df['num_star'] == num_stars]
    new_cores = new_coredf['with_star'].values.tolist()
    new_core = random.choice(new_cores)
    new_core = Chem.MolFromSmiles(new_core)
    combined = Chem.CombineMols(side_chains_combnind, new_core)
    ## 根据原子属性标记出接口虚拟原子
    isotope_dict = {1: []} # {虚拟原子Isotope:虚拟原子Idx},使用同位素虚拟原子来指定接口
    for atom in combined.GetAtoms():
        isotope = atom.GetIsotope()
        symbol = atom.GetSymbol()
        if isotope in isotope_dict and symbol == '*':
            isotope_dict[isotope].append(atom.GetIdx())
        if isotope not in isotope_dict and symbol == '*':
            isotope_dict[isotope] = [atom.GetIdx()]
    ## 相同标记的虚拟原子之间成键
    rw_mol = Chem.RWMol(combined)
    for key in isotope_dict.keys():
        isotope_key = isotope_dict[key]
        rw_mol.AddBond(isotope_key[0], isotope_key[1], Chem.BondType.SINGLE)
    modified_mol = rw_mol.GetMol() # 获取修改后的分子
    Chem.SanitizeMol(modified_mol) # 检查并清理分子
    ## 正则匹配式删除虚拟原子
    modified_smi = Chem.MolToSmiles(modified_mol)
    temp_file = (core, new_core, modified_mol)
    pattern = r'\[(\d+)\*\]|\(\[(\d+)\*\]\)'
    result_smi = re.sub(pattern, '', modified_smi)
    if Chem.MolFromSmiles(result_smi): # 检查替换后的分子是否有效，避免多价虚拟原子。
        result_mol = Chem.MolFromSmiles(result_smi)
        return result_smi, result_mol, temp_file
    else:
        result_mol = None
        return result_smi, result_mol, temp_file

def replacecore_many(smi: str, data_path: str, num: int) -> List[str]:
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    list_smi = []
    for _ in range(num):
        try:
            result_smi, result_mol, _ = replace_core_randomvision(mol, data_path)
            if result_mol is not None:
                list_smi.append(result_smi)
        except:
            pass
    return list_smi


if __name__ == "__main__":
    result_smi, result_mol = replace_core(
    'ClC1=CC(N2CCCC2)=C(C=CN3C4=CC=CC=C4)C3=C1', 
    'C12=CC=CC=C1NC=C2', 
    '[1*]C1=C([2*])N=C([3*])N=N1')
    print(result_smi)

