##### import business #####
import sympy as sp
import string, sys, itertools

##### run Python read in params #####
if len(sys.argv) > 1:
    print(f"Received input: spin-{sys.argv[1]}")
    spin_S = int(sys.argv[1])
else:
    print("insufficient input provided")


def gen_poly(spin_S):

    """
    Generates two polynomials:
        1. Poly A represents the sum of all coefficients of all possible diagrams 
           based on the input parameter spin_S under the no-offset assumption.
        2. Poly B represents the remainder polynomial that needs to be added to Poly A
           to ensure the correct overall polynomial is obtained when plus-offset is switched on.

    Parameters:
    spin_S (int): An integer representing the spin value.

    Returns: [for both Poly A and Poly B]
    sympy.Poly: A SymPy polynomial object representing the polynomial with coefficients 
                and powers derived from the constructed diagrams and their interactions.

    Description:
    1.  Construct a dictionary of coefficients for different power values based on spin_S.
        - `key_list` contains even numbers between 0 and 2*spin_S.
        - `value_list` contains coefficients, starting with 1 followed by symbolic coefficients.

    2.  Construct a dictionary (`mother_dict`) where keys are tuples representing all 
        possible combinations of values from `key_list` and values are the products of 
        corresponding coefficients from `coeff_dict`.

    3.  Constructs daughter diagrams from each mother diagram.
        - Each mother diagram generates up to four possible daughter diagrams by permuting 
          the positions of the labels.
        - `daughter_dict` contains these permutations as keys and their respective 
          coefficients (or negative coefficients) as values.

    4.  Constructs a dictionary (`polynomial_dict`) where:
        - Keys are powers of x derived from the sum of labels in the daughter diagrams.
        - Values are the accumulated coefficients for these powers.

    5.  Constructs the final polynomial expression using SymPy.
        - The polynomial is created by summing terms of the form x^power * coeff.

    6.  Based on Poly A, Poly B is constructed by extracting coefficients from Poly A involving more 
        than one algebraic coefficient in {a,b,c,d...}, e.g., a**2 or a*b but not a solo a, b etc.

    Example:
    >>> gen_polynomial(2)
    Poly A: Poly(x**4 - 2*a*x**3 + (a**2 + 2)*x**2 - 2*a*x + 1, x, domain='ZZ[a]')
    Poly B: 
    """

    spin_2S = spin_S * 2

    # create dict where key list is even numbers between 0 and spin_S
    # and value list is the coefficents 1, a, b, c, ...
    key_list = [spin_S*2 - i for i in range(0, spin_S+1, 2)]
    if len(key_list) > 2:
        value_list = list(sp.symbols(' '.join(string.ascii_lowercase[:len(key_list)-1])))
    else:
        value_list = [sp.symbols('a')]
    value_list.insert(0, 1)
    coeff_dict = dict(zip(key_list, value_list))

    # mother diagram dict: construct unique diagrams
    key_cartesian_product = list(itertools.product(key_list, key_list))
    mother_dict = {}
    for pair in key_cartesian_product:
        mother_dict[pair] = coeff_dict[pair[0]] * coeff_dict[pair[1]]

    # daughter diagram dict: each mother diagram generates {1,2,4} daughter diagrams
    # key are tuples (i,j,k,l): i,j are top labels y,x and k,l are bottom labels x,y
    daughter_dict = {}
    for pair, coeff in mother_dict.items():
        i,j,k,l = pair[0], spin_2S - pair[0], pair[1], spin_2S - pair[1]
        quad_permutes = list(set([(i,j,k,l), (j,i,k,l), (i,j,l,k), (j,i,l,k)]))        
        for quad_permute in quad_permutes:
            if (quad_permute[0] + quad_permute[2]) % 4 == 0:
                daughter_dict[quad_permute] = coeff
            else:
                daughter_dict[quad_permute] = -coeff

    # construct poly A dict: key is power of z, value is coeff
    # we transformed x = yd and z = d**2
    poly_A_dict = {}
    for quad_permute, coeff in daughter_dict.items():
        power = int((quad_permute[0] + quad_permute[2]) / 2)
        if power in poly_A_dict:
            poly_A_dict[power] += coeff
        else:
            poly_A_dict[power] = coeff
    poly_A_dict = dict(sorted(poly_A_dict.items(), reverse=True))

    # construct poly B dict
    poly_B_dict = gen_remainder_poly(value_list, poly_A_dict)

    # construct the polynomial in sympy: poly A and poly B
    poly_A, poly_B, x = 0, 0, sp.symbols('x')
    for power, coeff in poly_A_dict.items():
        poly_A += x**power * coeff

    for power, coeff in poly_B_dict.items():
        poly_B += x**power * coeff

    return sp.Poly(poly_A, x), sp.Poly(poly_B, x)

def gen_remainder_poly(value_list, poly_A_dict):
    
    remainder_dict = {}
    for key, coeff in poly_A_dict.items():
        if coeff == 1: continue
        for term in coeff.as_ordered_terms():
            if "-" in str(term):
                if -term not in value_list and term / -2 not in value_list:
                    if key in remainder_dict:
                        remainder_dict[key] += term
                    else:
                        remainder_dict[key] = term
            else:
                if term not in value_list and term / 2 not in value_list:
                    if key in remainder_dict:
                        remainder_dict[key] += term
                    else:
                        remainder_dict[key] = term


    return remainder_dict


def gen_poly_dict(polynomial):
    """
    Prints the terms and coefficients of a SymPy polynomial and returns a dictionary of these terms.

    Parameters:
    polynomial (sympy.Poly): A SymPy polynomial object to be printed and analyzed.

    Returns:
    dict: A dictionary where keys are the powers of the polynomial and values are the corresponding coefficients.

    Description:
    - Extracts terms and coefficients from the given polynomial.
    - Creates a dictionary with the power of the term as the key and the coefficient as the value.
    - Prints each power and its associated coefficient.
    - Returns the dictionary for further use.

    Example:
    >>> poly = Poly(x**4 - 2*a*x**3 + (a**2 + 2)*x**2 - 2*a*x + 1, x, domain='ZZ[a]')
    >>> print_polynomial(poly)
            4: 1
            3: -2*a
            2: a**2 + 2
            1: -2*a
            0: 1
            {4: 1, 3: -2*a, 2: a**2 + 2, 1: -2*a, 0: 1}
    """

    terms = polynomial.terms()
    coeff_dict = {term[0][0]: term[1] for term in terms}
    return coeff_dict

def factorize_poly(polynomial):
    """
    Factorizes a SymPy polynomial and prints the factors.

    Parameters:
    polynomial (sympy.Poly or sympy.Expr): A SymPy polynomial object or expression to be factorized.

    Returns:
    sympy.Expr: The factored polynomial expression.

    Description:
    - Uses SymPy's `factor` function to factorize the given polynomial.
    - Prints the factored form of the polynomial.
    - Returns the factored polynomial for further use.

    Example:
    >>> x = sp.symbols('x')
    >>> poly = sp.Poly(x**4 - 2*x**3 + x**2 - 2*x + 1, x)
    >>> factorized_poly = factorize_polynomial(poly)
    (x - 1)**2 * (x**2 + 1)
    """
    # Ensure the polynomial is in the form of a SymPy expression
    if isinstance(polynomial, sp.Poly):
        polynomial = polynomial.as_expr()
    
    # Factorize the polynomial
    factored_poly = sp.factor(polynomial)
    
    return factored_poly

# Generate the polynomials: poly A and poly B
poly_A, poly_B = gen_poly(spin_S)
coeff_dict_poly_A, coeff_dict_poly_B = gen_poly_dict(poly_A), gen_poly_dict(poly_B)

factor_poly_A, factor_poly_B = factorize_poly(poly_A), factorize_poly(poly_B)

def print_poly_dict(poly_dict):
    for key, value in coeff_dict_poly_A:
        print(f"{key}: {value}")

# Print the factored polynomial
sp.pprint(poly_A, use_unicode=True)
sp.pprint(factor_poly_A, use_unicode=True)
sp.pprint(poly_B, use_unicode=True)
sp.pprint(factor_poly_B, use_unicode=True)