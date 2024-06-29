#!/bin/python
# -*- coding: utf-8 -*-

def main():
    numbers = []
    with open("config.dat", "r") as input_file:
        for line in input_file:
            # 去除行首尾的空格和换行符
            line = line.strip()
            # 以#号为分隔符，获取#号前面的内容
            parts = line.split("#")[0].strip().split()
            for part in parts:
                try:
                    number = int(part)
                    numbers.append(str(number))
                except ValueError:
                    pass

    with open("config1.dat", "w") as output_file:
        output_file.write(" ".join(numbers))

if __name__ == "__main__":
    main()
