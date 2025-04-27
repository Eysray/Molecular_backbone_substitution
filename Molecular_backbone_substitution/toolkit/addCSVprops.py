import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import List, Dict, Tuple, Optional # 导入类型提示
# 定义一个函数，用于计算属性并添加到 CSV 文件中
def add_properties_to_csv(input_csv, output_csv):
    # 读取 CSV 文件
    df = pd.read_csv(input_csv)

    # 定义要计算的属性
    properties = {
        'LOGP': Descriptors.MolLogP, # 脂溶性
        # 'MW': Descriptors.MolWt,
        'TPSA': Descriptors.TPSA, # 拓扑极性表面积
        'HeavyAtomCount': Descriptors.HeavyAtomCount, # heavy atom
        'NumRotatableBonds': Descriptors.NumRotatableBonds, #旋转键
        'HBD': Descriptors.NumHDonors, #氢键供体
        'HBA': Descriptors.NumHAcceptors, #氢键受体
        'RingCount': Descriptors.RingCount, #环
        'FractionCSP3': Descriptors.FractionCSP3, #CSP3
        'NumValenceElectrons': Descriptors.NumValenceElectrons, #价电子
    }

    # 遍历每一行，计算属性
    for prop_name, prop_func in properties.items():
        values = []
        for smiles in df['SMILES']:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                value = prop_func(mol)
            else:
                value = None
            values.append(value)
        df[prop_name] = values

    # 保存为新的 CSV 文件
    df.to_csv(output_csv, index=False)
    
def add_smi_prop(smi: str) -> Optional[Dict]:
    # 读取 CSV 文件
    mol = Chem.MolFromSmiles(smi)
    # 定义要计算的属性
    properties_dict = {
        'LOGP': Descriptors.MolLogP, # 脂溶性
        'MW': Descriptors.MolWt,
        'TPSA': Descriptors.TPSA, # 拓扑极性表面积
        'HeavyAtomCount': Descriptors.HeavyAtomCount, # heavy atom
        'NumRotatableBonds': Descriptors.NumRotatableBonds, #旋转键
        'HBD': Descriptors.NumHDonors, #氢键供体
        'HBA': Descriptors.NumHAcceptors, #氢键受体
        'RingCount': Descriptors.RingCount, #环
        'FractionCSP3': Descriptors.FractionCSP3, #CSP3
        'NumValenceElectrons': Descriptors.NumValenceElectrons, #价电子
    }
    properties = {}
    # 遍历每一行，计算属性
    for prop_name, prop_func in properties_dict.items():
        resault = {
            'smiles': smi,
            'properties':properties,
        }
        if mol is not None:
            value = prop_func(mol)
            properties[prop_name] = value
    return resault


if __name__ == "__main__":
    input_csv = 'CL-amidine.csv'  # 替换为你的输入 CSV 文件
    output_csv = 'CL-amidine_addProps.csv'  # 替换为你的输出 CSV 文件
    add_properties_to_csv(input_csv, output_csv)
    print(add_smi_prop('C1=CC=CC=C1'))