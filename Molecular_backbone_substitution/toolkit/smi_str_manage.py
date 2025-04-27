import re, csv
from itertools import permutations, product
from typing import List # 导入类型提示
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm
'''本脚本用于处理SMILES字符串，包括：生成/清除虚拟原子；给有虚拟原子的SMILES字符串添加同位素编号；统计虚拟原子的数目；排列组合虚拟同位素的位置'''
'''
输入一个字符串，其中包含多个 [数字*] 部分，
例如 '[1*]C1=C([2*])N=C([3*])N=N1
输出所有可能的填充后的字符串，
'''
def fill_combinations(input_string: str) -> List[str]:
    '''生成排列组合的取代编号顺序'''
    # 使用正则表达式找到所有的 [数字*] 部分，并提取数字
    pattern = r'\[(\d+)\*\]'
    matches = re.findall(pattern, input_string)
    numbers = list(map(int, matches))  # 提取的数字列表
    # 生成所有可能的排列
    unique_perms = set(permutations(numbers))  # 使用集合去重（如果有重复数字）
    # 替换函数
    def replace_numbers(perm):
        # 将输入字符串拆分成部分，方便替换
        parts = re.split(r'(\[\d+\*\])', input_string)
        perm_iter = iter(perm)
        # 遍历 parts，遇到 [数字*] 时替换为当前排列的数字
        for i in range(len(parts)):
            if re.fullmatch(r'\[\d+\*\]', parts[i]):
                parts[i] = f'[{next(perm_iter)}*]'
        return ''.join(parts)
    # 生成所有填充后的字符串
    filled_strings = [replace_numbers(p) for p in unique_perms]
    return filled_strings


def remove_staratom(smi: str) -> str:
    #清除所有的虚拟原子
    pattern = r'\[(\d+)\*\]|\(\[(\d+)\*\]\)'
    smi = re.sub(pattern, '', smi)
    return smi

def from_txt_get_chemdraw_smicopy(txt_file_path: str) -> str:
    # 标准TXT文件 -> ChemDraw多复制字符串
    # 避免TXT空行
    with open(txt_file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        new_lines = []
        for line in lines:
            if line[-1] == '\n':  # 检查最后一个字符是否是换行符
                line = line[:-1] + '.' # 去除换行符
            new_lines.append(line)
    content = "".join(new_lines)  # 拼接成字符串
    cleaned = content.replace(" ", "") #清除空格
    return cleaned    


def process_chemdraw_smicopy(input_string, has_star = True) -> List[str]:
    # ChemDraw多复制字符串(如果含有不带同位素虚拟原子) -> 带同位素虚拟原子的SMILES字符串列表
    # 也接受分割好的含有[*]的SMILES列表
    if type(input_string) == str:
        smi_list = input_string.split('.') #按'.'分割为列表
    if type(input_string) == list:
        smi_list = input_string
    if has_star:
        add_num_list = []
        for smi in smi_list: #为每个虚拟原子'[*]'按顺序添加编号
            count = 1
            result = ""
            index = 0
            while index < len(smi):
                if smi[index:index + 3] == "[*]":
                    result += f"[{count}*]"
                    count += 1
                    index += 3
                else:
                    result += smi[index]
                    index += 1
            add_num_list.append(result)
        return add_num_list
    else:
        return smi_list

def count_stars(input_string: str) -> int:
    # 计算字符串中'*'的数量
    return input_string.count('*')

def gennerate_datacsv(str_list: List[str], csv_file_path: str, marker = None) -> None:
    # 带同位素虚拟原子的SMILES字符串列表 -> 带虚拟原子和不带虚拟原子的SMILES字符串CSV文件，以及虚拟原子数目
    if marker is not None:
        headers = ['SMILES', 'with_star', 'no_star', 'num_star', 'ring_count','logP', 'TPSA', 'HBD', 'HBA','marker']
    else:
        headers = ['SMILES', 'with_star', 'no_star', 'num_star', 'ring_count','logP', 'TPSA', 'HBD', 'HBA']
    with open(csv_file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for i in tqdm(range(len(str_list))):
            smi = str_list[i]
            real_smi = remove_staratom(smi)
            mol = Chem.MolFromSmiles(real_smi) #使用真实分子计算性质
            Chem.Kekulize(mol, clearAromaticFlags=True) # 清除芳香性,以匹配SMART匹配检索，注意：共振式无法检索
            if mol is not None:
                if marker is not None:
                    writer.writerow(
                        [smi, #带同位素虚拟原子的SMILES
                        smi, #带同位素虚拟原子的SMILES
                        real_smi, #不带同位素虚拟原子的SMILES，对应的完整SMILES
                        count_stars(smi), #虚拟原子数目
                        Descriptors.RingCount(mol), # 环数
                        Descriptors.MolLogP(mol), # 脂性
                        Descriptors.TPSA(mol), # 拓扑表面积
                        Descriptors.NumHDonors(mol), # 氢键供体
                        Descriptors.NumHAcceptors(mol), # 氢键受体
                        marker[i], # 标记
                        ])
                else:
                    writer.writerow(
                        [smi, #带同位素虚拟原子的SMILES
                        smi, #带同位素虚拟原子的SMILES
                        real_smi, #不带同位素虚拟原子的SMILES
                        count_stars(smi), #虚拟原子数目
                        Descriptors.RingCount(mol), # 环数
                        Descriptors.MolLogP(mol), # 脂性
                        Descriptors.TPSA(mol), # 拓扑表面积
                        Descriptors.NumHDonors(mol), # 氢键供体
                        Descriptors.NumHAcceptors(mol), # 氢键受体
                        ])
            else:
                print(f'errorsmi:{smi}')

def generate_replacements(original_str: str, target='[H]', replacement='[*]') -> List[str]:
    # 原始SMILES字符串（显示全部氢原子） -> 排列组合的带虚拟原子的SMILES字符串列表（包括同位重复）
    #随即替换某个字符串部分为另一个字符串部分
    # Find all the positions where target occurs
    positions = []
    start = 0
    target_len = len(target)
    while True:
        pos = original_str.find(target, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    n = len(positions)
    replacements = []
    # Generate all possible combinations of replacing or not replacing each target
    for bits in tqdm(product([False, True], repeat=n)):
        temp_str = list(original_str)
        # Iterate from the end to the beginning to avoid messing up positions
        for i in range(n-1, -1, -1):
            if bits[i]:
                pos = positions[i]
                temp_str[pos:pos+target_len] = replacement
        replacements.append(''.join(temp_str))
    return replacements

def add_star_atom(smi: str) -> List[str]:
    # 虚拟原子衍生化
    #输入一个不带虚拟原子/氢原子的SMILES字符串，输出带虚拟原子不带氢原子的SMILES字符串列表
    mol = Chem.MolFromSmiles(smi)
    Chem.Kekulize(mol, clearAromaticFlags=True) # 清除芳香性,以适应杂环衍生
    mol_h = Chem.AddHs(mol)  # 添加氢原子
    smi_h = Chem.MolToSmiles(mol_h)  # 转换回SMILES
    star_h_list = generate_replacements(smi_h)
    star_list = []
    pattern = r'\[H\]|\(\[H\]\)'
    for star_h in tqdm(star_h_list): # 清除所有的氢原子
        star = re.sub(pattern, '', star_h)
        star_list.append(star)
    return star_list

if __name__ == "__main__":
    # smi = '[1*]C1=C([2*])N=C([3*])N=N1'  # 输入字符串
    # filled_strings = fill_combinations(smi)  # 生成所有填充后的字符串
    # removeed = remove_staratom(smi)
    # print(filled_strings)  # 输出所有填充后的字符串
    # print(removeed)
    # input_string = 'CN[*].CN(C)[*].C[N+](C)(C)[*].[*]C1=C([*])C=CC=C1.[*]C2=C([*])C=C([*])C=C2'
    # print(process_chemdraw_smicopy(input_string))
    x = add_star_atom('CC1=CC=CS1')
    mols = [Chem.MolFromSmiles(smi) for smi in x]
    import rdkit.Chem.Draw as Draw
    print(x)
    img = Draw.MolsToGridImage(mols, molsPerRow=10)
    img.save('test.png')
