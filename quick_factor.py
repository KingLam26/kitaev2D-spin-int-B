import sys
import sympy as sp

# Define all variables
z = sp.symbols('z')
a1, a2, a3 = sp.symbols('a1 a2 a3')
b1, b2, b3 = sp.symbols('b1 b2 b3')
c1, c2, c3 = sp.symbols('c1 c2 c3')

# Define palindromal polynomial coefficients
poly_coeff_dict = {4: [1, -2*a1, 2+a1**2], 
				   6: [1, -2*a1, 2*a2 + a1**2, -(2*a1+2*a1*a2), 2+2*a1**2+a2**2],
				   8: [1, -2*a1, 2*a2 + a1**2, -(2*a3+2*a1*a2), 2*a2+2*a1*a3+a2**2, -(2*a1+2*a1*a2+2*a2*a3), 2+2*a1**2+2*a2**2+a3**2],
				   5: [1, -(2*a1+2), 3+a1**2-2*a1],
				   7: [1, -(2*a1-2), 3-4*a1+2*a2+a1**2, -(-4+6*a1-2*a2-2*a1**2+2*a1*a2), 5-6*a1+2*a2+3*a1**2-2*a1*a2+a2**2],
				   #1.5: [a2**2, -(2*a2 + 2*a1*a2**2), 1+a1**2*a2**2+2*a1*a2+2*a2**2],
				   1.5: [a1**2, -(2*a1 + 2*a1*a1**2), 1+a1**2*a1**2+2*a1*a1+2*a1**2],
				   2.5: [a2**2, -(2*a1*a2 + 2*a1*a2**2), a1**2 + 2*a2**3 + 2*a2+4*a1**2*a2+a1**2*a2**2, -(2*a1 + 4*a1*a2 + 2*a1**3 + 4*a1*a2**2 + 2*a1**3*a2 + 2*a2**3*a1), 1+ 4*a1**2 + 4*a2**2 + 4*a1**2*a2 + 4*a1**2*a2**2 + a1**4 + a2**4]
				   }				   
# Generate palindromal polynomial
def gen_PP(poly_coeff_list, spin):
	poly = 0
	count = len(poly_coeff_list)
	deg = (count-1) * 2
	for i in range(count):
		if i*2 != deg:
			poly += poly_coeff_list[i] * z**i + poly_coeff_list[i] * z**(deg-i)
		else:
			poly += poly_coeff_list[i] * z**i
	return poly

# Factor palindromal polynomial
if len(sys.argv) > 1:
	spin = float(sys.argv[1]) if '.' in sys.argv[1] else int(sys.argv[1])

poly = gen_PP(poly_coeff_dict[spin], spin)
poly_factor = sp.factor(poly)
print(poly_factor)