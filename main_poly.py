from utils_poly import *

diff_print = False
poly_print_A, poly_print_B = True, True
quo_A_print, quo_B_print = True, True
factor_A_print, factor_B_print = True, True


##### run Python read in params #####
if len(sys.argv) > 1:
    print(f"Received input: spin-{sys.argv[1]}")
    spin_S = int(sys.argv[1])
else:
    print("insufficient input provided")


##### main analysis #####
# Generate the polynomials: poly A and poly B
poly_A, poly_B = gen_poly(spin_S)

if poly_print_A:
    print(f"\nPoly A: \n")
    print_poly_dict(gen_poly_dict(poly_A))
    
if poly_print_B:
    print(f"\nPoly B: \n")
    print_poly_dict(gen_poly_dict(poly_B))

# Differentiate Poly A
first_diff_poly_A = diff_poly(poly_A)
second_diff_poly_A = diff_poly(first_diff_poly_A)

first_diff_dict_poly_A = gen_poly_dict(first_diff_poly_A)
second_diff_dict_poly_A = gen_poly_dict(second_diff_poly_A)

if diff_print:
    print(f"\nPoly A first diff: \n")
    print_poly_dict(first_diff_dict_poly_A)
    print(f"\nPoly A second diff: \n")
    print_poly_dict(second_diff_dict_poly_A)

# Divide Poly A by (x-1)**2 if S is odd, else no change
x = sp.symbols('x')
if spin_S % 2 != 0:
    denom_A = (x-1)**2
else:
    denom_A = 1
quo_poly_A, remain_poly_A = sp.div(poly_A.as_expr(), denom_A, domain = "QQ")
quo_poly_A = sp.Poly(quo_poly_A, x)
if remain_poly_A != 0:
    print("remainder!")
    quit()

if quo_A_print:
    print(f"\nPoly A quotient: \n")
    print_poly_dict(gen_poly_dict(quo_poly_A))

# Divide Poly B by x**2 always, and (x-1)**2 if S odd
if spin_S % 2 != 0:
    denom_B = x**2 * (x-1)**2
else:
    denom_B = x**2
quo_poly_B, remain_poly_B = sp.div(poly_B.as_expr(), denom_B, domain = "QQ")
quo_poly_B = sp.Poly(quo_poly_B, x)
if remain_poly_B != 0:
    print("remainder!")
    quit()

if quo_B_print:
    print(f"\nPoly B quotient: \n")
    print_poly_dict(gen_poly_dict(quo_poly_B))

# factorize Poly A
factors_A_list = sp.factor_list(quo_poly_A)
if factor_A_print:
    print(f"\nPoly A factor: \n")
    print_poly_dict(gen_poly_dict(factors_A_list[1][0][0]))

# factorize Poly B
factors_B_list = sp.factor_list(quo_poly_B)
if factor_B_print:
    print(f"\nPoly B factor: \n")
    print_poly_dict(gen_poly_dict(factors_B_list[1][0][0]))

quit()