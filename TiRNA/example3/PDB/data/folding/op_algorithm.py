import os
import subprocess

def main():
    # 读取配置文件
    with open("config1.dat", "r") as fp:
        lines = fp.readlines()
        # 取最后一行（原C代码的行为）
        last_line = lines[-1].strip()
        ca, cb, cc, cd, ce, cf = map(int, last_line.split())
    
    # 切换到 folding 目录
    
    # 根据 ca 的值执行不同的脚本
    if ca == 1:
        subprocess.run(["bash", "run_remc.sh"])
    else:
        subprocess.run(["bash", "run_sa.sh"])

if __name__ == "__main__":
    main()
