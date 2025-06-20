""" import sys
import gzip

with gzip.open(sys.argv[1],"r") as f:
	for line in f:
		lsplit = line.rstrip().split("\t")
		if lsplit[7]:
			line = lsplit[0] + "\t" + lsplit[7]
			print(line)
 """

import sys
import gzip

print(">>> 开始读取文件 <<<")
with gzip.open(sys.argv[1], "rt") as f:  # 使用文本模式读取

       
    for line in f:
        lsplit = line.rstrip().split("\t")
        if len(lsplit) > 7 and lsplit[7]:  # 检查索引存在性
            new_line = f"{lsplit[0]}\t{lsplit[7]}"  # 使用f-string更清晰
            print(new_line)

""" import sys
import gzip

with gzip.open(sys.argv[1], "rt") as f:
    line_counter = 0  # 行计数器
    print(">>> 开始读取文件 <<<")
    
    for line in f:
        line_counter += 1
        print(f"\n=== 正在处理第 {line_counter} 行 ===")
        print("原始内容:", repr(line))  # 显示换行符等特殊字符

        # 处理字段分割
        cleaned_line = line.rstrip()
        lsplit = cleaned_line.split("\t")
        print(f"字段分割结果（共 {len(lsplit)} 列）: {lsplit}")

        # 检查字段有效性
        if len(lsplit) > 7:
            print(f"第8字段存在，内容为: {repr(lsplit[7])}")
            if lsplit[7]:
                new_line = f"{lsplit[0]}\t{lsplit[7]}"
                print(f"★ 有效行 - 重组内容: {new_line}")
                print(new_line)  # 原打印逻辑
            else:
                print("× 第8字段为空值，跳过")
        else:
            print(f"! 警告：本行只有 {len(lsplit)} 列，不满足条件")

    print(f"\n>>> 处理完成，共扫描 {line_counter} 行 <<<") """