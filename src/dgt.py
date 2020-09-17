from gaussian import GaussianInteger
from params import p, g, g_inv

def dgt_gentleman_sande(A):
    n = len(A)

    m = n//2
    while m >= 1:
        for j in range(m):
            a = pow(g, (j * n)//(2 * m), p)
            i = j
            while i < n:
                xi = A[i]
                xim = A[i + m]
                
                A[i] = (xi + xim) % p
                A[i + m] = (a * (xi - xim)) % p

                i = i + 2*m
        m = m//2

    return A

def idgt_cooley_tukey(A):
    n = len(A)

    m = 1
    while m < n:
        for j in range(m):
            a = pow(g_inv, int((j * n)/(2 * m)), p)

            i = j
            while i < n:
                xi = A[i]
                xim = A[i + m]

                A[i] = (xi + a * xim) % p
                A[i + m] = (xi - a * xim) % p

                i = i + 2*m
        m = 2*m 

    return A