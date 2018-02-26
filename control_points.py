#!/usr/bin/env python
# Code from https://www.particleincell.com/wp-content/uploads/2012/06/bezier-spline.js

def print_points(name, points):
    points = map(lambda p: "({:.5f}pt,{:.5f}pt)".format(p[0],p[1]), points)
    print("{}: [{}]".format(name, ", ".join(points)))

points = [
    (0,1),
    (1,0),
    (1.5,3),
    (2,1),
    (3,1)
]

K = list(map(lambda x: (x[0]*28.45274, x[1]*28.45274), points))
#K.reverse()
print("points: {}".format(points))
print_points("K", K)
n = len(K) - 1

a = [0] * n
b = [0] * n
c = [0] * n
r = [(0,0)] * n
p = [(0,0)] * n
q = [(0,0)] * n

# Left most segment
a[0] = 0
b[0] = 2
c[0] = 1
r[0] = (K[0][0]+2*K[1][0],
        K[0][1]+2*K[1][1])

# Internal segments
for i in range(1, n-1):
    a[i] = 1
    b[i] = 4
    c[i] = 1
    r[i] = (4*K[i][0]+2*K[i+1][0],
            4*K[i][1]+2*K[i+1][1])

# Right segment
a[n-1] = 2
b[n-1] = 7
c[n-1] = 0
r[n-1] = (8*K[n-1][0]+K[n][0],
          8*K[n-1][1]+K[n][1])

# Solve Ax=b via Thomas's algorithm
for i in range (1, n):
    m = a[i]/b[i-1]
    print("a[i]/b[i-1] = m: {}/{}={}".format(a[i], b[i-1], m))
    b[i] = b[i] - m*c[i-1]
    r[i] = (r[i][0] - m*r[i-1][0],
            r[i][1] - m*r[i-1][1])

print("a: {}\nb: {}\nc: {}".format(a,b,c))
print_points("r", r)

# Compute the first control points p
p[n-1] = (r[n-1][0]/b[n-1],
          r[n-1][1]/b[n-1])
for i in range(n-2, -1, -1):
    p[i] = ((r[i][0] - c[i]*p[i+1][0])/b[i],
            (r[i][1] - c[i]*p[i+1][1])/b[i])

# Compute the second control points q
for i in range(0, n-1):
    q[i] = (2*K[i+1][0] - p[i+1][0],
            2*K[i+1][1] - p[i+1][1])
q[n-1] = ((K[n][0] + p[n-1][0])/2,
          (K[n][1] + p[n-1][1])/2)

for i in range(0, n):
    print(".. controls ({:.5f}pt,{:.5f}pt) and ({:.5f}pt,{:.5f}pt) .. ({:.5f}pt,{:.5f}pt)" \
            .format(p[i][0], p[i][1],
                    q[i][0], q[i][1],
                    K[i+1][0], K[i+1][1]))
