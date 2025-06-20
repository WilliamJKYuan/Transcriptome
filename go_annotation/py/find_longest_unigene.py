import sys
import re
x = sys.argv[1]
def readfasta(filename):
    fa = open(filename,'r')
    res = {}
    ID = ''
    for line in fa:
        if line.startswith('>'):
            ID = line.strip('\n')
            res[ID] = ''
        else:
            res[ID] +=line.strip('\n')
    return res
res = readfasta(x)
uniq = {}
regex = re.compile(r'=\d+')
longest = 0
for k,v in res.items():
    i = k.split(' ',2)
    title = (i[0]).split('_')
    title = '_'.join(title[:4]) + '\n'
    v = [v]
    if title not in uniq:
        uniq[title] = v
    else:
        uniq[title] += v
max_seq = {}
for k,v in uniq.items():
    seq = max(v, key = len)
    max_seq[k] =seq
rlt = str(max_seq)
rlt = rlt.replace("{'", '')
rlt = rlt.replace("\\n': '", "\n")
rlt = rlt.replace("', '", "\n")
rlt = rlt.replace("'}","\n").rstrip()
print(rlt)