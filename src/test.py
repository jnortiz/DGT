import unittest
import random
from gaussian import GaussianInteger
from dgt import idgt_cooley_tukey, dgt_gentleman_sande
from params import N, p, kth_root_of_i, inv_kth_root_of_i, invkmodp

NRUNS = 10**2

fold   = lambda a: [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
unfold = lambda a: [a[j].re for j in range(N//2)] + [a[j].imag for j in range(N//2)]

def gen_polynomial_modp(length):
    x = []

    for _ in range(length):
        x.append(random.randrange(0,p))
    
    return x

def schoolbook_mul(a, b):
    assert len(a) == len(b)
    
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i] * b[j] * (-1)**(int((i+j)//float(N)))
            c[(i + j) % N] = (c[(i + j) % N] + v) % p
    
    return c

def poly_mul_gs_then_ct(x, y):

    # Folding the input signal into N//2 Gaussian integers
    x_folded = fold(x)

    # Applying the right-angle convolution
    x_folded = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

    # Folding the input signal into N//2 Gaussian integers
    y_folded = fold(y)

    # Applying the right-angle convolution
    y_folded = [y_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

    # Computing the transform on both operands
    x_dgt = dgt_gentleman_sande(x_folded)
    y_dgt = dgt_gentleman_sande(y_folded)

    # Point-wise multiplication in Z_p[i]
    xy_dgt = [(xi * yi) % p for xi, yi in zip(x_dgt, y_dgt)]

    # Computing the inverse DGT transform
    xy = idgt_cooley_tukey(xy_dgt)

    # Removing the twisting factors and scaling by k^-1
    xy = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(xy)]

    # Unfolding the output signal
    xy = unfold(xy)

    return xy

class TestDGT(unittest.TestCase):

    def test_dgt_gs_then_ct(self):
        print("\nDGT transformation (Gentleman-Sande then Cooley-Tukey)")

        for _ in range(NRUNS):
            
            # Generating at random an n-degree polynomial with coefficients modulo p
            x = gen_polynomial_modp(N)

            # Folding the input signal into N//2 Gaussian integers
            x_folded = fold(x)

            # Applying the right-angle convolution
            x_folded_twisted = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

            # Computing both forward and inverse transforms, in this order
            y_folded_twisted = idgt_cooley_tukey(dgt_gentleman_sande(x_folded_twisted))

            # Removing the twisting factors and scaling by k^-1
            y_folded = [(xi * invkmodp * inv_kth_root_of_i[i]) % p for i, xi in enumerate(y_folded_twisted)]

            # Unfolding the output signal
            y = unfold(y_folded)
           
            self.assertEqual(x, y)

    def test_mul_gs_then_ct(self):
        print("\nPolynomial multiplication using DGT (Gentleman-Sande then Cooley-Tukey)")

        for _ in range(NRUNS):
            
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)

            ab = schoolbook_mul(a, b)
            c = poly_mul_gs_then_ct(a, b)

            self.assertEqual(c, ab)

    def test_twisting(self):
        print("\nRight-angle convolution (folding and twisting procedures)")

        for _ in range(NRUNS):

            # Generating at random an n-degree polynomial with coefficients modulo p
            x = gen_polynomial_modp(N)

            # Folding the input signal into N//2 Gaussian integers
            x_folded = fold(x)

            # Applying the right-angle convolution
            x_folded_twisted = [x_folded[i] * kth_root_of_i[i] % p for i in range(N//2)]

            # Removing the twisting factors and scaling by k^-1
            y_folded = [(xi * inv_kth_root_of_i[i]) % p for i, xi in enumerate(x_folded_twisted)]

            # Unfolding the output signal
            y = unfold(y_folded)
            
            self.assertEqual(x, y)

if __name__ == '__main__':
    unittest.main()