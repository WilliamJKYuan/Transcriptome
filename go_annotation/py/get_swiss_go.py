import sys 
import pandas as pd
import progress

print("Start")

print("start importing file 1")
id_df = pd.read_csv(sys.argv[1],sep='\t',header=None)
print("success import file 1")


print("start importing file 2")
trinity_df = pd.read_csv(sys.argv[2],sep='\t',header=None) 
print("success import file 2")

id_df.columns = ['swiss_id','go_id']
trinity_df.columns = ['trinity_id', 'swiss_id']

df = pd.merge(id_df, trinity_df, on = 'swiss_id', how = 'inner')

df2 = df[['trinity_id','go_id']]
df2.to_csv('Results/genes.go.annotation',index=None, sep = '\t')
print("output success")
