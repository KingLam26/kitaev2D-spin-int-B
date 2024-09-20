import sympy as sp

# Define the variables
x, z, a, b, c = sp.symbols('x z a b c')

# Define the polynomial
poly_4 = x**4 - 2*a*x**3 + (2+a**2)*x**2 -2*a*x + 1

poly_6 = z**8 - 2*a*z**7 + (2*b + a**2)*z**6 - (2*a + 2*a*b)*z**5 + (2 + 2*a**2 + b**2)*z**4 - (2*a + 2*a*b)*z**3 + (2*b + a**2)*z**2 - 2*a*z + 1

poly_8 = z**12 - 2*a*z**11 + (2*b+a**2)*z**10 - (2*c+2*a*b)*z**9 + (2*b+b**2+2*a*c)*z**8 - (2*a+2*a*b+2*b*c)*z**7 + (2+2*a**2+2*b**2+c**2)*z**6
poly_8 += 1 - 2*a*z**1 + (2*b+a**2)*z**2 - (2*c+2*a*b)*z**3 + (2*b+b**2+2*a*c)*z**4 - (2*a+2*a*b+2*b*c)*z**5

# Factorize the polynomial
poly_4_factor = sp.factor(poly_4)
poly_6_factor = sp.factor(poly_6)
poly_8_factor = sp.factor(poly_8)

# Display the factorized form
print(poly_4_factor)
print(poly_6_factor)
print(poly_8_factor)
