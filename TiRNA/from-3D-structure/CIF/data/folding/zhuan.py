import os

def cif_to_ch_dat():
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
    
    input_cif = "initial.cif"
    output_file = "ch.dat"
    
    try:
        with open(input_cif, 'r') as f_in:
            lines = f_in.readlines()
        
        with open(output_file, 'w') as f_out:
            atom_count = 0
            residue_data = {}  # 按残基存储原子数据
            second_residue_p_coord = None  # 记录第二个残基的P原子坐标
            
            # 解析CIF文件的原子坐标部分
            in_atom_section = False
            atom_lines = []
            column_indices = {}  # 存储各字段的列索引
            
            # 找到原子坐标部分并确定字段位置
            for i, line in enumerate(lines):
                line = line.strip()
                
                # 查找原子坐标部分的开始
                if line.startswith('_atom_site.'):
                    in_atom_section = True
                    # 记录字段位置
                    field_name = line.replace('_atom_site.', '').split()[0]
                    column_indices[field_name] = len(column_indices)
                    continue
                
                # 收集原子坐标行
                if in_atom_section and line and not line.startswith('_') and not line.startswith('#') and not line.startswith('loop_'):
                    atom_lines.append(line)
                elif in_atom_section and (line.startswith('_') or line.startswith('loop_')):
                    # 遇到新的部分，停止收集
                    if atom_lines:  # 如果已经收集到原子数据，就停止
                        break
            
            # 解析原子坐标
            for line in atom_lines:
                parts = line.split()
                if len(parts) < 6:
                    continue
                
                try:
                    # 根据字段索引获取相应的值
                    # 原子名称可能有不同的字段名
                    atom_name = ""
                    for field in ['label_atom_id', 'auth_atom_id', 'type_symbol']:
                        if field in column_indices:
                            atom_name = parts[column_indices[field]]
                            break
                    
                    if not atom_name:
                        continue
                    
                    # 残基名称
                    res_name = ""
                    for field in ['label_comp_id', 'auth_comp_id']:
                        if field in column_indices:
                            res_name = parts[column_indices[field]]
                            break
                    
                    # 链ID
                    chain_id = "A"
                    for field in ['label_asym_id', 'auth_asym_id']:
                        if field in column_indices:
                            chain_id = parts[column_indices[field]]
                            break
                    
                    # 残基序号
                    res_seq = 1
                    for field in ['label_seq_id', 'auth_seq_id']:
                        if field in column_indices:
                            res_seq_str = parts[column_indices[field]]
                            # 处理可能带引号的残基序号
                            res_seq_str = res_seq_str.strip("'\"")
                            if res_seq_str and res_seq_str != '.':
                                try:
                                    res_seq = int(float(res_seq_str))  # 处理带小数的序号
                                except:
                                    res_seq = 1
                            break
                    
                    # 坐标字段
                    x = float(parts[column_indices.get('Cartn_x', len(parts)-3)])
                    y = float(parts[column_indices.get('Cartn_y', len(parts)-2)])
                    z = float(parts[column_indices.get('Cartn_z', len(parts)-1)])
                    
                    # 只处理核苷酸
                    valid_residues = ['A', 'G', 'C', 'U', 'T', 'DA', 'DG', 'DC', 'DT']
                    if res_name not in valid_residues:
                        continue
                    
                    # 标准化残基名称（处理DNA残基）
                    if res_name in ['DA', 'DG', 'DC', 'DT']:
                        base_name = res_name[1]  # 去掉D前缀
                    else:
                        base_name = res_name
                    
                    # 创建残基键
                    residue_key = f"{res_seq}_{base_name}_{chain_id}"
                    
                    # 初始化残基数据
                    if residue_key not in residue_data:
                        residue_data[residue_key] = {
                            'res_name': base_name,
                            'res_seq': res_seq,
                            'chain_id': chain_id,
                            'atoms': {}
                        }
                    
                    # 收集需要的原子：P、C4'、N1/N9
                    # 处理可能的原子名称变体
                    if atom_name in ['P']:  # 磷酸基团原子
                        residue_data[residue_key]['atoms']['P'] = (x, y, z)
                        # 记录第二个残基的P原子坐标
                        if len(residue_data) == 2 and second_residue_p_coord is None:
                            second_residue_p_coord = (x, y, z)
                    
                    # C4'原子的各种可能名称
                    elif atom_name in ["C4'", "C4*", "C4"]:  # 糖环C4'原子
                        residue_data[residue_key]['atoms']["C4'"] = (x, y, z)
                    
                    elif atom_name in ['N1', 'N9']:  # 碱基特征原子
                        residue_data[residue_key]['atoms'][atom_name] = (x, y, z)
                        
                except (IndexError, ValueError, KeyError):
                    continue
            
            if not residue_data:
                return
            
            # 按残基序号排序并输出
            sorted_residues = sorted(residue_data.values(), key=lambda x: (x['chain_id'], x['res_seq']))
            
            # 检查第一个残基是否有P原子
            first_residue_has_p = 'P' in sorted_residues[0]['atoms'] if sorted_residues else False
            
            # 如果第一个残基没有P原子，且第二个残基有P原子，则创建一个虚拟P原子
            if sorted_residues and not first_residue_has_p and second_residue_p_coord:
                x, y, z = second_residue_p_coord
                virtual_x = x + 2.0
                virtual_y = y + 2.0
                virtual_z = z + 2.0
                first_res_seq = sorted_residues[0]['res_seq']
                
                f_out.write(f"{atom_count+1} {first_res_seq} P {virtual_x:.6f} {virtual_y:.6f} {virtual_z:.6f} {radius_map['P']:.6f} {charge_map['P']:.6f} 0.000000\n")
                atom_count += 1
            
            # 正常输出所有残基的3个珠子
            for residue in sorted_residues:
                res_name = residue['res_name']
                res_seq = residue['res_seq']
                atoms = residue['atoms']
                
                # 输出P原子（如果第一个残基没有P原子且这是第一个残基，则跳过，因为已经创建了虚拟P原子）
                if 'P' in atoms and not (residue == sorted_residues[0] and not first_residue_has_p):
                    x, y, z = atoms['P']
                    f_out.write(f"{atom_count+1} {res_seq} P {x:.6f} {y:.6f} {z:.6f} {radius_map['P']:.6f} {charge_map['P']:.6f} 0.000000\n")
                    atom_count += 1
                
                # 输出C4'原子（在输出文件中标记为S）
                if "C4'" in atoms:
                    x, y, z = atoms["C4'"]
                    f_out.write(f"{atom_count+1} {res_seq} S {x:.6f} {y:.6f} {z:.6f} {radius_map['S']:.6f} {charge_map['S']:.6f} 0.000000\n")
                    atom_count += 1
                else:
                    # 如果没有C4'原子，尝试使用其他糖环原子作为替代
                    for atom_name in atoms:
                        if atom_name.startswith('C') and atom_name not in ['P', 'N1', 'N9']:
                            x, y, z = atoms[atom_name]
                            f_out.write(f"{atom_count+1} {res_seq} S {x:.6f} {y:.6f} {z:.6f} {radius_map['S']:.6f} {charge_map['S']:.6f} 0.000000\n")
                            atom_count += 1
                            break
                
                # 输出碱基特征原子（N1/N9）
                base_atom_name = base_atoms.get(res_name)
                if base_atom_name and base_atom_name in atoms:
                    x, y, z = atoms[base_atom_name]
                    f_out.write(f"{atom_count+1} {res_seq} {res_name} {x:.6f} {y:.6f} {z:.6f} {radius_map.get(res_name, 2.150):.6f} {charge_map.get(res_name, 0.000):.6f} 0.000000\n")
                    atom_count += 1
                elif res_name in base_atoms:
                    # 如果没有特征原子，尝试使用其他碱基原子作为替代
                    for atom_name in atoms:
                        if atom_name.startswith('N') or atom_name.startswith('C'):
                            x, y, z = atoms[atom_name]
                            f_out.write(f"{atom_count+1} {res_seq} {res_name} {x:.6f} {y:.6f} {z:.6f} {radius_map.get(res_name, 2.150):.6f} {charge_map.get(res_name, 0.000):.6f} 0.000000\n")
                            atom_count += 1
                            break
            
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
                atom_count += 1
        
    except FileNotFoundError:
        pass
    except Exception:
        pass

def auto_convert():
    """自动检测文件类型并进行转换"""
    if os.path.exists("initial.pdb"):
        # 这里可以调用原来的PDB转换函数
        pass
    elif os.path.exists("initial.cif"):
        cif_to_ch_dat()
    else:
        pass

if __name__ == "__main__":
    auto_convert()
