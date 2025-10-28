import os
import glob
from Bio.PDB import PDBParser, MMCIFIO
import warnings
warnings.filterwarnings("ignore")

def convert_pdb_to_cif_in_current_folder():
    """将当前文件夹中的所有 PDB 文件转换为 CIF 格式"""
    
    # 获取当前文件夹中的所有 PDB 文件
    pdb_files = glob.glob("*.pdb")
    
    if not pdb_files:
        return
    
    parser = PDBParser(QUIET=True)
    io = MMCIFIO()
    
    for pdb_file in pdb_files:
        try:
            # 获取文件名（不含扩展名）
            filename = os.path.splitext(pdb_file)[0]
            output_file = f"{filename}.cif"
            
            # 读取并转换
            structure = parser.get_structure(filename, pdb_file)
            io.set_structure(structure)
            io.save(output_file)
            
        except Exception:
            # 静默失败，不输出任何信息
            pass

# 直接运行
convert_pdb_to_cif_in_current_folder()
