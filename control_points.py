#!/usr/bin/env python
# Code from https://www.particleincell.com/wp-content/uploads/2012/06/bezier-spline.js

def print_points(name, points):
    points = map(lambda p: "({:.5f}pt,{:.5f}pt)".format(p[0],p[1]), points)
    print("{}: [{}]".format(name, ", ".join(points)))

def spline_through(points):
    K = list(map(lambda x: (x[0]*28.45274, x[1]*28.45274), points))
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
    
    print(r"\draw[red, ->] ({:.5f}pt,{:.5f}pt)".format(K[0][0], K[0][1]))
    for i in range(n):
        print(".. controls ({:.5f}pt,{:.5f}pt) and ({:.5f}pt,{:.5f}pt) .. ({:.5f}pt,{:.5f}pt)" \
                .format(p[i][0], p[i][1],
                        q[i][0], q[i][1],
                        K[i+1][0], K[i+1][1]))
    print(";")
    for i in range(n):
        print(r"\fill[green] ({:.5f}pt,{:.5f}pt) circle (1pt);".format(p[i][0], p[i][1]))
        print(r"\fill[red]   ({:.5f}pt,{:.5f}pt) circle (1pt);".format(q[i][0], q[i][1]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[i][0], K[i][1], p[i][0], p[i][1]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[i+1][0], K[i+1][1], q[i][0], q[i][1]))

# a_i = 1;
# b_i = 4
# c_i = 1
# b_{n-1} = 3
# alpha = beta = gamma = 1
# Solve Ax = rx, Ay = ry
def cyclic2(rx, ry):
    n = len(rx)
    assert n == len(ry)

    x = [-1e9] * n
    y = [-1e9] * n
    z = [-1e9] * n
    bb = [-1e9] * n

    # Since alpha = beta = gamma = 1,
    # b[0] = b[n-1] = 3 and u[0] = u[n-1] = 1
    bb[0] = 3.0
    x[0] = rx[0]
    y[0] = ry[0]
    z[0] = 1
    for i in range(1, n-1):
        bb[i] = 4.0 - 1.0/bb[i-1]       # bb[i] =  b[i] - c[i-1]*a[i]/bb[i-1]
        x[i] = rx[i] - x[i-1]/bb[i-1]   #  x[i] = rx[i] - x[i-1]*a[i]/bb[i-1]
        y[i] = ry[i] - y[i-1]/bb[i-1]   #  y[i] = ry[i] - y[i-1]*a[i]/bb[i-1]
        z[i] = 0.0 - z[i-1]/bb[i-1]     #  z[i] =  u[i] - z[i-1]*a[i]/bb[i-1]
    bb[n-1] = 3.0 - 1/bb[n-2]           # bb[n-1] =  b[n-1] - c[n-2]*a[n-1]/bb[n-2]
    x[n-1] = rx[n-1] - x[n-2]/bb[n-2]   #  x[n-1] = rx[n-1] - c[n-2]*a[n-1]/bb[n-2]
    y[n-1] = ry[n-1] - y[n-2]/bb[n-2]   #  y[n-1] = ry[n-1] - c[n-2]*a[n-1]/bb[n-2]
    z[n-1] = 1.0 - z[n-2]/bb[n-2]

    x[n-1] = x[n-1]/bb[n-1]             # x[n-1] = x[n-1]/bb[n-1]
    y[n-1] = y[n-1]/bb[n-1]             # y[n-1] = y[n-1]/bb[n-1]
    z[n-1] = z[n-1]/bb[n-1]             # z[n-1] = z[n-1]/bb[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (x[i] - x[i+1])/bb[i]    # x[i] = (x[i] - c[i]*x[i+1])/bb[i]
        y[i] = (y[i] - y[i+1])/bb[i]    # y[i] = (y[i] - c[i]*y[i+1])/bb[i]
        z[i] = (z[i] - z[i+1])/bb[i]    # z[i] = (z[i] - c[i]*z[i+1])/bb[i]

    d = 1 + z[0] + z[n-1]
    fx = (x[0] + x[n-1])/d
    fy = (y[0] + y[n-1])/d

    for i in range(n):
        x[i] -= fx*z[i]
        y[i] -= fy*z[i]
    return x, y

def closed_spline_through3(points):
    K = list(map(lambda x: (x[0]*28.45274, x[1]*28.45274), points))
    n = len(points)

    rx = [0] * n
    ry = [0] * n
    Qx = [0] * n
    Qy = [0] * n

    for i in range(n):
        j = (i+1) % n
        rx[i] = 4*K[i][0] + 2*K[j][0]
        ry[i] = 4*K[i][1] + 2*K[j][1]

    Px, Py = cyclic2(rx, ry)
    for i in range(n):
        j = (i+1) % n
        Qx[i] = 2*K[j][0] - Px[j]
        Qy[i] = 2*K[j][1] - Py[j]

    print(r"\draw[blue, ->] ({:.5f}pt,{:.5f}pt)".format(K[0][0], K[0][1]))
    for i in range(n):
        j = (i+1) % n
        print(".. controls ({:.5f}pt,{:.5f}pt) and ({:.5f}pt,{:.5f}pt) .. ({:.5f}pt,{:.5f}pt)" \
                .format(Px[i], Py[i],
                        Qx[i], Qy[i],
                        K[j][0], K[j][1]))
    print(";")
    for i in range(n):
        j = (i+1) % n
        print(r"\fill[green] ({:.5f}pt,{:.5f}pt) circle (1pt);".format(Px[i], Py[i]))
        print(r"\fill[red]   ({:.5f}pt,{:.5f}pt) circle (1pt);".format(Qx[i], Qy[i]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[i][0], K[i][1], Px[i], Py[i]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[j][0], K[j][1], Qx[i], Qy[i]))


# Solve Ax = r where A is tridiagonal and defined by bands a (below the main
# diagonal), b (the main diagonal), and c (above the main diagonal).
# a[0] and c[n-1] are not used.
def tridiag(a, b, c, r):
    n = len(a)
    assert n == len(b)
    assert n == len(c)
    assert n == len(r)

    x = [0.0] * n
    gam = [-1e9] * n

    bet = b[0]
    x[0] = r[0]/bet
    for i in range(1,n):
        gam[i] = c[i-1]/bet
        bet = b[i] - a[i]*gam[i]
        x[i] = (r[i] - a[i]*x[i-1])/bet

    for i in range(n-2, -1, -1):
        x[i] -=  gam[i+1]*x[i+1]
    return x

def cyclic(a, b, c, alpha, beta, r, gamma=1.0):
    n = len(a)
    assert n == len(b)
    assert n == len(c)
    assert n == len(r)
    
    bb = b.copy()
    bb[0] = b[0] - gamma
    bb[n-1] = b[n-1] - alpha*beta/gamma
    x = tridiag(a, bb, c, r)

    u = [0] * n
    u[0] = gamma
    u[n-1] = alpha
    z = tridiag(a, bb, c, u)

    factor = (x[0] + beta*x[n-1]/gamma)/(1 + z[0] + beta*z[n-1]/gamma)

    for i in range(n):
        x[i] -= factor*z[i]
    return x

def closed_spline_through2(points):
    K = list(map(lambda x: (x[0]*28.45274, x[1]*28.45274), points))

    n = len(points)
    a = [1.0] * n
    b = [4.0] * n
    c = [1.0] * n
    rx = [0.0] * n
    ry = [0.0] * n

    for i in range(n):
        j = (i+1) % n
        rx[i] = 4*K[i][0] + 2*K[j][0]
        ry[i] = 4*K[i][1] + 2*K[j][1]
    x = cyclic(a, b, c, 1, 1, rx)
    y = cyclic(a, b, c, 1, 1, ry)

    P = list(zip(x,y))
    Q = [(0,0)] * n
    for i in range(n):
        j = (i+1) % n
        Q[i] = (2*K[j][0] - P[j][0],
                2*K[j][1] - P[j][1])
    print(r"\draw[blue, ->] ({:.5f}pt,{:.5f}pt)".format(K[0][0], K[0][1]))
    for i in range(n):
        j = (i+1) % n
        print(".. controls ({:.5f}pt,{:.5f}pt) and ({:.5f}pt,{:.5f}pt) .. ({:.5f}pt,{:.5f}pt)" \
                .format(P[i][0], P[i][1],
                        Q[i][0], Q[i][1],
                        K[j][0], K[j][1]))
    print(";")
    for i in range(n):
        j = (i+1) % n
        print(r"\fill[green] ({:.5f}pt,{:.5f}pt) circle (1pt);".format(P[i][0], P[i][1]))
        print(r"\fill[red]   ({:.5f}pt,{:.5f}pt) circle (1pt);".format(Q[i][0], Q[i][1]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[i][0], K[i][1], P[i][0], P[i][1]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[j][0], K[j][1], Q[i][0], Q[i][1]))


def closed_spline_through(points):
    K = list(map(lambda x: (x[0]*28.45274, x[1]*28.45274), points))
    print_points("K", K)

    n = len(K)
    a = [1] * n
    b = [4] * n
    c = [1] * n
    rx = [0] * n
    ry = [0] * n
    x = [0] * n
    y = [0] * n
    z = [0] * n
    u = [0] * n
    
    for i in range(n-1):
        rx[i] = 4*K[i][0] + 2*K[i+1][0]
        ry[i] = 4*K[i][1] + 2*K[i+1][1]
    rx[n-1] = 4*K[n-1][0] + 2*K[0][0]
    ry[n-1] = 4*K[n-1][1] + 2*K[0][1]

    alpha = a[0]
    beta = c[n-1]
    gamma = 1
    b[0] = b[0] - gamma
    b[n-1] = b[n-1] - alpha*beta/gamma

    u[0] = gamma
    u[n-1] = alpha

    # Solve Ax = rx, Ay = ry and Az=u via Thomas's algorithm
    x[0] = rx[0] / b[0]
    y[0] = ry[0] / b[0]
    for i in range(1,n):
        m = a[i]/b[i-1]
        b[i] = b[i] - m*c[i-1]
        rx[i] = rx[i] - m*rx[i-1]
        ry[i] = ry[i] - m*ry[i-1]
        u[i] = u[i] - m*u[i-1]
    
    # Back substitution
    x[n-1] = rx[n-1]/b[n-1]
    y[n-1] = ry[n-1]/b[n-1]
    z[n-1] = u[n-1]/b[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (rx[i] - c[i]*x[i+1])/b[i]
        y[i] = (ry[i] - c[i]*y[i+1])/b[i]
        z[i] = (u[i] - c[i]*z[i+1])/b[i]

    # Compute the multiplicative factor for x and y
    fx = (x[0] + beta*x[n-1]/gamma)/(1 + z[0] + beta*z[n-1]/gamma)
    fy = (y[0] + beta*y[n-1]/gamma)/(1 + z[0] + beta*z[n-1]/gamma)

    p = [(0,0)] * n
    q = [(0,0)] * n
    for i in range(n):
        j = i % n
        p[i] = (x[i] - fx*z[i],
                y[i] - fy*z[i])
        q[i] = (2*K[j][0] - p[j][0],
                2*K[j][1] - p[j][1])

    print(r"\draw[blue, ->] ({:.5f}pt,{:.5f}pt)".format(K[0][0], K[0][1]))
    for i in range(n-1):
        print(".. controls ({:.5f}pt,{:.5f}pt) and ({:.5f}pt,{:.5f}pt) .. ({:.5f}pt,{:.5f}pt)" \
                .format(p[i][0], p[i][1],
                        q[i][0], q[i][1],
                        K[i+1][0], K[i+1][1]))
    print(".. controls ({:.5f}pt,{:.5f}pt) and ({:.5f}pt,{:.5f}pt) .. ({:.5f}pt,{:.5f}pt);" \
            .format(p[n-1][0], p[n-1][1],
                    q[n-1][0], q[n-1][1],
                    K[0][0], K[0][1]))
    for i in range(n):
        if i == n-1:
            j = 0
        else:
            j = i+1
        print(r"\fill[green] ({:.5f}pt,{:.5f}pt) circle (1pt);".format(p[i][0], p[i][1]))
        print(r"\fill[red]   ({:.5f}pt,{:.5f}pt) circle (1pt);".format(q[i][0], q[i][1]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[i][0], K[i][1], p[i][0], p[i][1]))
        print(r"\draw ({:.5f}pt,{:.5f}pt) -- ({:.5f}pt,{:.5f}pt);".format(K[j][0], K[j][1], q[i][0], q[i][1]))

if __name__ == '__main__':
    points = [
        (0,1),
        (3,1),
        (2,3),
        (1.5,0),
    ]
    #print("Spline through points:")
    #spline_through(points)
    #print("Closed spline through points:")
    closed_spline_through2(points)
    closed_spline_through3(points)

# vim: set sw=4 sts=4 ts=8 et tw=0:
