from typing import List, Dict, Tuple, Optional # 导入类型提示
from rdkit import Chem # 导入环境
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG # 导入绘图工具

def draw_mol_svg(mol: Chem.Mol, figursize: Tuple, filename = None, AtomIndices = False, BondIndices = False) -> None:
    """绘制分子结构输出SVG文件"""
    d2d = rdMolDraw2D.MolDraw2DSVG(*figursize)
    if AtomIndices: # 显示原子索引
        d2d.drawOptions().addAtomIndices = True
    if BondIndices: # 显示键索引
        d2d.drawOptions().addBondIndices = True
    d2d.DrawMolecule(mol) # 绘制分子
    d2d.FinishDrawing() # 必须调用FinishDrawing()方法，否则不会显示
    svg_fig = SVG(d2d.GetDrawingText()) # SVG对象在Jupyter Notebook中可以直接显示（不指定'svg_fig ='的情况下）
    if filename: # 保存到文件
        with open(filename, 'w') as f:
            f.write(d2d.GetDrawingText())
    else: # 直接显示
        return svg_fig