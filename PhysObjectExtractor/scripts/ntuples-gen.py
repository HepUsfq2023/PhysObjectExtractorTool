import json

f = open("nevts_recid.txt", "r")
lines = f.readlines()
f.close()
size=len(lines)

variables=[]
for i in range(0, size):
    aux=lines[i].split()
    variables.append([aux[0], int(aux[1]), aux[2] ])

def findRecid(recid):
    result=" "
    for j in range(0,size):
        if(recid==variables[j][0]):
            result= variables[j]    
    return result      

data = {}
data['data'] =  {}
data['ttbar'] = {}
data['wjets'] =  {}

#24132
data['data']['nominal'] = {}
data['data']['nominal']['files']=[]
data['data']['nominal']['files'].append({
    'path': findRecid("24132")[2],
    'nevts': findRecid("24132")[1]
})



#19359
data['ttbar']['nominal']={}
data['ttbar']['nominal']['files']=[]
data['ttbar']['nominal']['files'].append({
    'path': findRecid("19359")[2],
    'nevts': findRecid("19359")[1]
})
data['ttbar']['nominal']['nevts_total']=findRecid("19359")[1]

#19949
data['wjets']['ME_var']={}
data['wjets']['ME_var']['files']=[]
data['wjets']['ME_var']['files'].append({
    'path': findRecid("19949")[2],
    'nevts': findRecid("19949")[1]
})
data['wjets']['ME_var']['nevts_total']=findRecid("19949")[1]

with open('ntuples.json', 'w') as files:
    json.dump(data, files)
