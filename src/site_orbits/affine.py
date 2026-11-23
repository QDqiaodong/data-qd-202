from fractions import Fraction

def determinant(m):
    return sum(m[0][i]*(m[1][(i+1)%3]*m[2][(i+2)%3]-m[1][(i+2)%3]*m[2][(i+1)%3]) for i in range(3))

def key(matrix,translation):return (tuple(tuple(row) for row in matrix),tuple(t%1 for t in translation))

def compose(left,right):
    a,u=left;b,v=right
    matrix=tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    t=tuple((sum(a[i][j]*v[j] for j in range(3))+u[i])%1 for i in range(3))
    return matrix,t

def apply(op,point):
    matrix,t=op
    return tuple((sum(row[j]*point[j] for j in range(3))+t[i])%1 for i,row in enumerate(matrix))

IDENTITY=key(((1,0,0),(0,1,0),(0,0,1)),(Fraction(0),)*3)
