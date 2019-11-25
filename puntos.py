ny=5
nx=2
b=4
c=1
yn=[]

for i in range(0,nx) :
    for j in range(0,ny) :
        yn.append([i*c+0.25*c,b*j/(ny-1)-(b/2)])
print(yn)
