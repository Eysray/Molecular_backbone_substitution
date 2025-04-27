import toolkit as tk
import pandas as pd
# 所有分子要清除芳香性
def main1():
    #单取代库生成
    # txt文档最后一行不要有空格，不要有换行，避免多余的'.'添加
    txt_file_path = r'data\chemdraw.txt'
    csv_file_path = r'data\smi_data.csv'
    cleaned = tk.from_txt_get_chemdraw_smicopy(txt_file_path) # 标准TXT文件 -> ChemDraw多复制字符串
    smi_withstar = tk.process_chemdraw_smicopy(cleaned) # ChemDraw多复制字符串(如果含有不带同位素虚拟原子) -> 带同位素虚拟原子的SMILES字符串列表
    tk.gennerate_datacsv(smi_withstar, csv_file_path) # 带同位素虚拟原子的SMILES字符串列表 -> 带虚拟原子和不带虚拟原子的SMILES字符串CSV文件，以及虚拟原子数目

    
def main2():
    # 环取代库生成
    # txt文档最后一行不要有空格，不要有换行，避免多余的'.'添加
    txt_file_path = r'data\ring.txt'
    csv_file_path = r'data\ring_data.csv'
    cleaned = tk.from_txt_get_chemdraw_smicopy(txt_file_path) # 标准TXT文件 -> ChemDraw多复制字符串
    smi_withstar = tk.process_chemdraw_smicopy(cleaned, has_star = False) # 无虚拟原子，分割多复制字符串
    star_list = []
    for smi in smi_withstar:
        tmplist = tk.add_star_atom(smi) #会产生不含有虚拟原子的结构
        tmplist = tk.process_chemdraw_smicopy(tmplist) #虚拟原子编码
        star_list += tmplist # 虚拟原子衍生化
    tk.gennerate_datacsv(star_list, csv_file_path) # 带同位素虚拟原子的SMILES字符串列表 -> 带虚拟原子和不带虚拟原子的SMILES字符串CSV文件，以及虚拟原子数目     
    
def main3():
    # txt文档最后一行不要有空格，不要有换行，避免多余的'.'添加
    smi = 'C1=CC=CN1'
    csv_file_path = r'data\ring_debug_data.csv'
    tmplist = tk.add_star_atom(smi) #会产生不含有虚拟原子的结构
    tmplist = tk.process_chemdraw_smicopy(tmplist) #虚拟原子编码
    tk.gennerate_datacsv(tmplist, csv_file_path) # 带同位素虚拟原子的SMILES字符串列表 -> 带虚拟原子和不带虚拟原子的SMILES字符串CSV文件，以及虚拟原子数目   

def main4():
    # 读取两个 CSV 文件
    df1 = pd.read_csv(r'data\ring_data.csv')
    df2 = pd.read_csv(r'data\smi_data.csv')
    # 按行合并（上下堆叠）
    merged_df = pd.concat([df1, df2], ignore_index=True)
    # 保存结果
    merged_df.to_csv(r'data\merged.csv', index=False)  
        
if __name__ == "__main__":
    main1()
    main2()
    main3()
    main4()