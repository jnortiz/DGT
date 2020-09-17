from gaussian import GaussianInteger

modinv = lambda y,p:pow(y,p-2,p)

p = 0xFFFFFFFF00000001
N = 512
k = N//2

kthroots = {
	16: GaussianInteger(18446462319560032257, 18014402800254976),
	32: GaussianInteger(18446181093691752453, 18446181110871097349),
	64: GaussianInteger(1512311202781865130, 475568752833222537),
	128: GaussianInteger(16359014788679420321, 9463965103391723449),
	256: GaussianInteger(1100507988529617178,13061373484646814047),
	512: GaussianInteger(14193207829130683548, 15237665595666387946),
	1024: GaussianInteger(11605794709690543498, 8056921223206052788),
	2048: GaussianInteger(1372650973816680927, 13948452854352226361),
	4096: GaussianInteger(13224239242624218619, 468760267013759104),
	8192: GaussianInteger(1345080250393555844, 4827528560175661949)
}

invkthroots = {
	16: GaussianInteger(72057594037944320, 4611686018427388160),
	32: GaussianInteger(1152912708648048640, 1152912708379604992),
	64: GaussianInteger(13162941959907671889, 12615130990074419764),
	128: GaussianInteger(2065520003433100932, 16643184082144809425),
	256: GaussianInteger(1012656084342873654, 2372108221600941182),
	512: GaussianInteger(16889751086295376116, 18147277524140977891),
	1024: GaussianInteger(5650167340200895162, 12297659996166169359),
	2048: GaussianInteger(15611101961941206979, 13033551770418110946),
	4096: GaussianInteger(9806235711004847349, 3105654936161609),
	8192: GaussianInteger(11434313317047606333, 11715719925817632455)
}

# Primitive root of p
PROOTS = {
	0xFFFFFFFF00000001:7
}

# Computing the k-th primitive root of p
r = PROOTS[p] # Primitive root of p
assert (p-1) % k == 0
n = (p-1)//k

g = int(pow(r, n, p))
g_inv = modinv(g, p)

assert pow(g, k, p) == 1 # k-th primitive root of unity
assert pow(g_inv, k, p) == 1 # Inverse of the k-th primitive root of unity

# Computing the powers of the k-th root of i mod p
assert pow(kthroots[N//2], N//2) % p == GaussianInteger(0, 1)
assert (kthroots[N//2] * invkthroots[N//2]) % p == GaussianInteger(1, 0)

kth_root_of_i = [pow(kthroots[N//2], i) % p for i in range(N//2)]
inv_kth_root_of_i = [pow(invkthroots[N//2], i) % p for i in range(N//2)]

# Scalar k^-1 mod p used in to scale the output signal in backward DGT
invkmodp = modinv(k, p)