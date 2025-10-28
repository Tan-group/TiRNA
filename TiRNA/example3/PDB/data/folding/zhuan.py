import os

def pdb_to_ch_dat():
    # 定义残基类型到半径和电荷的映射
    radius_map = {
        'P': 2.110,
        'S': 1.800,  # 对应C4'
        'G': 2.150,
        'C': 2.150,
        'A': 2.150,
        'U': 2.150
    }
    
    charge_map = {
        'P': -1.000,
        'S': 0.000,
        'G': 0.000,
        'C': 0.000,
        'A': 0.000,
        'U': 0.000
    }
    
    # 碱基对应的特征原子
    base_atoms = {
        'A': 'N9',  # 嘌呤
        'G': 'N9',  # 嘌呤
        'C': 'N1',  # 嘧啶
        'U': 'N1',  # 嘧啶
        'T': 'N1'   # 胸腺嘧啶（DNA）
    }
    
    input_pdb = "initial.pdb"
    output_file = "ch.dat"
    
    try:
        with open(input_pdb, 'r') as f_in:
            lines = f_in.readlines()
        
        with open(output_file, 'w') as f_out:
            atom_count = 0
            residue_data = {}  # 按残基存储原子数据
            second_residue_p_coord = None  # 记录第二个残基的P原子坐标
            
            # 解析PDB文件的原子坐标部分
            for line in lines:
                line = line.strip()
                
                # 只处理ATOM记录
                if not line.startswith('ATOM'):
                    continue
                
                try:
                    # 解析PDB格式的固定字段
                    atom_serial = int(line[6:11].strip())  # 原子序号
                    atom_name = line[12:16].strip()  # 原子名称
                    res_name = line[17:20].strip()  # 残基名称
                    chain_id = line[21:22].strip()  # 链ID
                    res_seq = int(line[22:26].strip())  # 残基序号
                    x = float(line[30:38].strip())  # X坐标
                    y = float(line[38:46].strip())  # Y坐标
                    z = float(line[46:54].strip())  # Z坐标
                    
                    # 只处理核苷酸
                    if res_name not in ['A', 'G', 'C', 'U', 'T', 'DA', 'DG', 'DC', 'DT']:
                        continue
                    
                    # 标准化残基名称（处理DNA残基）
                    if res_name in ['DA', 'DG', 'DC', 'DT']:
                        base_name = res_name[1]  # 去掉D前缀
                    else:
                        base_name = res_name
                    
                    # 创建残基键
                    residue_key = f"{res_seq}_{base_name}"
                    
                    # 初始化残基数据
                    if residue_key not in residue_data:
                        residue_data[residue_key] = {
                            'res_name': base_name,
                            'res_seq': res_seq,
                            'atoms': {}
                        }
                    
                    # 收集需要的原子：P、C4'、N1/N9
                    if atom_name == 'P':
                        residue_data[residue_key]['atoms']['P'] = (x, y, z)
                        # 记录第二个残基的P原子坐标
                        if len(residue_data) == 2 and second_residue_p_coord is None:
                            second_residue_p_coord = (x, y, z)
                    elif atom_name == "C4'":
                        residue_data[residue_key]['atoms']["C4'"] = (x, y, z)
                    elif atom_name in ['N1', 'N9']:
                        residue_data[residue_key]['atoms'][atom_name] = (x, y, z)
                        
                except (IndexError, ValueError) as e:
                    # print(f"解析PDB行时出错: {line}")
                    # print(f"错误: {e}")
                    continue
            
            # print(f"找到 {len(residue_data)} 个核苷酸残基")
            
            # 按残基序号排序并输出
            sorted_residues = sorted(residue_data.values(), key=lambda x: x['res_seq'])
            
            # 检查第一个残基是否有P原子
            first_residue_has_p = 'P' in sorted_residues[0]['atoms'] if sorted_residues else False
            
            # 如果第一个残基没有P原子，且第二个残基有P原子，则创建一个虚拟P原子
            if sorted_residues and not first_residue_has_p and second_residue_p_coord:
                # 用第二个残基的P原子坐标加2.0创建第一个P原子
                x, y, z = second_residue_p_coord
                virtual_x = x + 2.0
                virtual_y = y + 2.0
                virtual_z = z + 2.0
                first_res_seq = sorted_residues[0]['res_seq']
                
                f_out.write(f"{atom_count+1} {first_res_seq} P {virtual_x:.6f} {virtual_y:.6f} {virtual_z:.6f} {radius_map['P']:.6f} {charge_map['P']:.6f} 0.000000\n")
                atom_count += 1
                # print(f"创建虚拟P原子用于第一个残基，坐标: ({virtual_x:.6f}, {virtual_y:.6f}, {virtual_z:.6f})")
            
            # 正常输出所有残基的3个珠子
            for residue in sorted_residues:
                res_name = residue['res_name']
                res_seq = residue['res_seq']
                atoms = residue['atoms']
                
                # print(f"处理残基 {res_seq} {res_name}: 找到原子 {list(atoms.keys())}")
                
                # 输出P原子（如果第一个残基没有P原子且这是第一个残基，则跳过，因为已经创建了虚拟P原子）
                if 'P' in atoms and not (res_seq == sorted_residues[0]['res_seq'] and not first_residue_has_p):
                    x, y, z = atoms['P']
                    f_out.write(f"{atom_count+1} {res_seq} P {x:.6f} {y:.6f} {z:.6f} {radius_map['P']:.6f} {charge_map['P']:.6f} 0.000000\n")
                    atom_count += 1
                    # print(f"  写入P原子: ({x:.6f}, {y:.6f}, {z:.6f})")
                
                # 输出C4'原子（在输出文件中标记为S）
                if "C4'" in atoms:
                    x, y, z = atoms["C4'"]
                    f_out.write(f"{atom_count+1} {res_seq} S {x:.6f} {y:.6f} {z:.6f} {radius_map['S']:.6f} {charge_map['S']:.6f} 0.000000\n")
                    atom_count += 1
                    # print(f"  写入C4'原子: ({x:.6f}, {y:.6f}, {z:.6f})")
                
                # 输出碱基特征原子（N1/N9）
                base_atom_name = base_atoms.get(res_name)
                if base_atom_name and base_atom_name in atoms:
                    x, y, z = atoms[base_atom_name]
                    f_out.write(f"{atom_count+1} {res_seq} {res_name} {x:.6f} {y:.6f} {z:.6f} {radius_map.get(res_name, 2.150):.6f} {charge_map.get(res_name, 0.000):.6f} 0.000000\n")
                    atom_count += 1
                    # print(f"  写入{base_atom_name}原子: ({x:.6f}, {y:.6f}, {z:.6f})")
                elif res_name in base_atoms:
                    # print(f"  警告: 残基 {res_name} 未找到特征原子 {base_atoms[res_name]}")
                    pass
            
            # 在最后添加一个额外的P原子
            if sorted_residues:
                last_residue = sorted_residues[-1]
                if 'P' in last_residue['atoms']:
                    x, y, z = last_residue['atoms']['P']
                elif second_residue_p_coord:
                    x, y, z = second_residue_p_coord
                else:
                    x, y, z = 0.0, 0.0, 0.0
                
                new_x = x + 2.0
                new_y = y + 2.0
                new_z = z + 2.0
                new_res_seq = last_residue['res_seq'] + 1
                f_out.write(f"{atom_count+1} {new_res_seq} P {new_x:.6f} {new_y:.6f} {new_z:.6f} {radius_map['P']:.6f} {charge_map['P']:.6f} 0.000000\n")
                # print(f"已添加额外的P原子，坐标: ({new_x:.6f}, {new_y:.6f}, {new_z:.6f})")
                atom_count += 1
            
            # print(f"总共写入 {atom_count} 个原子")
        
        # print(f"PDB粗粒化转换完成！输出文件: {output_file}")
        
    except FileNotFoundError:
        # print(f"错误：找不到输入文件 {input_pdb}")
        # print("请确保 initial.pdb 文件存在于当前目录中")
        pass
    except Exception as e:
        # print(f"转换过程中出现错误: {e}")
        # import traceback
        # traceback.print_exc()
        pass

def auto_convert():
    """自动检测文件类型并进行转换"""
    if os.path.exists("initial.pdb"):
        # print("检测到 initial.pdb 文件，进行PDB转换...")
        pdb_to_ch_dat()
    elif os.path.exists("initial.cif"):
        # print("检测到 initial.cif 文件，进行MMCIF转换...")
        # 这里可以调用原来的MMCIF转换函数
        # print("MMCIF转换功能暂未在此版本实现")
        pass
    else:
        # print("错误：未找到 initial.pdb 或 initial.cif 文件")
        # print("请确保输入文件存在于当前目录中")
        pass

if __name__ == "__main__":
    auto_convert()
