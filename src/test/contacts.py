from numpy import *
import linx

gro = linx.gro('bhp.gro')
h = [i for i,j in enumerate(gro[0]) if j['name'][0]!='H']
rs = []
lr = -1
for i,a in enumerate(gro[0]):
    if i in h:
        r = a['resid']
        if r!=lr: rs.append([]) 
        rs[-1].append(i)
        lr = r

from numpy.linalg import norm
def md(x1,x2):
    m = 9999
    for i in x1:
        for j in x2:
            m = min(m,norm(i-j))

    return m

x = gro[1]
ds = []
N = len(rs)
for i in range(N):
    for j in range(i+2,N):
        ds.append(md(x[rs[i]], x[rs[j]]))

        
savetxt('cts.txt',ds)
