import sys

class GaussianInteger(dict):
    re = 0
    imag = 0
    p = None

    def __init__(self, re, imag = 0, p = None):
        assert type(re) == int
        assert type(imag) == int

        self.re = re
        self.imag = imag

        if type(p) == int:
            self.p = p
        elif type(p) == float:
            raise Exception(
                "Weird p! Received " + str(p)
                )
        dict.__init__(self, re = re, imag = imag)

    def conjugate(self):
        return GaussianInteger(self.re, -self.imag)

    def __valid_operand(self, b):
        if type(b) == int:
            b = GaussianInteger(b, 0, 1)
        elif type(b) == float:
            raise Exception(
                "Error! We shouldn't be dealing with floats. Found " + str(b)
                )
        elif not isinstance(b, GaussianInteger):
            raise Exception(
                "Don't know how to deal with %s" % b
                )
        return b

    def __repr__(self):
        return "GaussianInteger(%d, %d)" % (self.re, self.imag)

    def __eq__(self, b):
        b = self.__valid_operand(b)

        return self.re == b.re and self.imag == b.imag

    def __ne__(self, b):
        return not(self == b
            )

    def __add__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re + b.re),
            (self.imag + b.imag),
            self.p
            )

    def __radd__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re + b.re),
            (self.imag + b.imag),
            self.p
            )

    def __sub__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re - b.re),
            (self.imag - b.imag),
            self.p
            )

    def __rsub__(self, b):
        b = self.__valid_operand(b)

        return GaussianInteger(
            (self.re - b.re),
            (self.imag - b.imag),
            self.p
            )

    def __mul__(self, b):
        b = self.__valid_operand(b)
        # https://stackoverflow.com/questions/19621686/complex-numbers-product-using-only-three-multiplications
        # 
        # S1=ac,S2=bd, and S3=(a+b)(c+d). Now you can compute the results as 
        # A=S1−S2 and B=S3−S1−S2.
        # 
        s1 = self.re * b.re
        s2 = self.imag * b.imag
        s3 = (self.re + self.imag) * (b.re + b.imag) 

        return GaussianInteger(
            (s1 - s2),
            (s3 - s1 - s2),
            self.p 
            )

    def __rmul__(self, b):
        b = self.__valid_operand(b)

        s1 = self.re * b.re
        s2 = self.imag * b.imag
        s3 = (self.re + self.imag) * (b.re + b.imag) 

        return GaussianInteger(
            (s1 - s2),
            (s3 - s1 - s2),
            self.p
            )

    def __pow__(self, b):
        if(b == 0):
            return GaussianInteger(1, 0, self.p)

        exp = bin(b)[3:]
        value = self

        for i in range(len(exp)):
            if self.p is not None:
                value = value * value % self.p
                if(exp[i:i+1] == '1'):
                    value = value * self % self.p
            else:
                value = value * value
                if(exp[i:i+1] == '1'):
                    value = value * self
        return value

    # About divisions and remainders:
    # https://math.stackexchange.com/questions/889809/calculating-the-reminder-when-dividing-complex-numbers
    def __floordiv__(self, b):
        if type(b) == int:
            assert b != 0

            return GaussianInteger(
                self.re // b,
                self.imag // b,
                self.p
                )
        else:
            assert isinstance(b, GaussianInteger)

            q = complex(self.re, self.imag) / complex(b.re, b.imag)
            return GaussianInteger(int(round(q.real)), int(round(q.imag)), self.p)

    def __mod__(self, b):
        if type(b) == int:
            assert b != 0

            return GaussianInteger(
                int(self.re % b),
                int(self.imag % b),
                self.p
                )
        else:
            assert isinstance(b, GaussianInteger)
            return self - (self // b) * b
            
    def norm(self):
        return self.re*self.re + self.imag*self.imag

    def xgcd(self,other):
        quot=GaussianInteger(0,0,self.p); a1=GaussianInteger(1,0,self.p); b1=GaussianInteger(0,0,self.p); a2=GaussianInteger(0,0,self.p)
        b2=GaussianInteger(1,0,self.p); a = self; b = other

        if(b.norm() > a.norm()):
            a,b = b,a  # Swap a and b - need to start with a>b
            a1,b1,a2,b2 = a2,b2,a1,b1 # Swap (a1,b1) with (a2,b2)
        while (1):
            quot = a // b
            a %= b
            a1 -= quot*a2; b1 -= quot*b2
            if (a == 0):
                return b, a2, b2
            quot = b // a
            b %= a
            a2 -= quot*a1; b2 -= quot*b1
            if (b == GaussianInteger(0,0,self.p)):
                return a, a1, b1

    def gcd(self,other):
        return self.xgcd(other)[0]