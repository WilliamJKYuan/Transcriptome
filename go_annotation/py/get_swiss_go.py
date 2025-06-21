import sys
import gzip

print(">>> 开始读取文件 <<<")
with gzip.open(sys.argv[1], "rt") as f:  # 使用文本模式读取

       
    for line in f:
        lsplit = line.rstrip().split("\t")
        if len(lsplit) > 7 and lsplit[7]:  # 检查索引存在性
            new_line = f"{lsplit[0]}\t{lsplit[7]}"  # 使用f-string更清晰
            print(new_line)
