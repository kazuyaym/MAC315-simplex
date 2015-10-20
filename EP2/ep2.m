A = [1,-2,1,0
1,-2,0,1]

m = rows (A)
n = columns (A)

c = [-1,-3,0,0]'
c = [0,0,1,1]'
b = [10, 20]'
x = [0,0,10,20]'

[ind v] = simplex(A,b,c,m,n,x);
