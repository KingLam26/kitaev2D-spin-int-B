##### import business #####
import sympy as sp
import string, sys, itertools

##### run Python read in params #####
if len(sys.argv) > 1:
    print(f"Received input: spin-{sys.argv[1]}")
    spin_S = int(sys.argv[1])
else:
    print("insufficient input provided")


def gen_polynomial(spin_S):

    """
    Generates a polynomial representing the sum of all coefficients of all possible diagrams
    based on the input parameter spin_S.

    Parameters:
    spin_S (int): An integer representing the spin value. This value is used to construct
                  the effective polynomial and diagrams.

    Returns:
    sympy.Poly: A SymPy polynomial object representing the polynomial with coefficients 
                and powers derived from the constructed diagrams and their interactions.

    Description:
    1. Constructs a dictionary of coefficients for different power values based on spin_S.
       - `key_list` contains even numbers between 0 and 2*spin_S.
       - `value_list` contains coefficients, starting with 1 followed by symbolic coefficients.

    2. Creates a dictionary (`mother_dict`) where keys are tuples representing all 
       possible combinations of values from `key_list` and values are the products of 
       corresponding coefficients from `coeff_dict`.

    3. Generates daughter diagrams from each mother diagram.
       - Each mother diagram generates up to four possible daughter diagrams by permuting 
         the positions of the labels.
       - `daughter_dict` contains these permutations as keys and their respective 
         coefficients (or negative coefficients) as values.

    4. Constructs a dictionary (`polynomial_dict`) where:
       - Keys are powers of x derived from the sum of labels in the daughter diagrams.
       - Values are the accumulated coefficients for these powers.

    5. Constructs the final polynomial expression using SymPy.
       - The polynomial is created by summing terms of the form x^power * coeff.

    Example:
    >>> gen_polynomial(2)
    Poly(x**4 - 2*a*x**3 + (a**2 + 2)*x**2 - 2*a*x + 1, x, domain='ZZ[a]')
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

    # polynomial dict: key is power of z, value is coeff
    # we transform x = yd and z = d**2
    polynomial_dict = {}
    for quad_permute, coeff in daughter_dict.items():
        power = int((quad_permute[0] + quad_permute[2]) / 2)
        if power in polynomial_dict:
            polynomial_dict[power] += coeff
        else:
            polynomial_dict[power] = coeff
    polynomial_dict = dict(sorted(polynomial_dict.items(), reverse=True))

    # construct the polynomial in sympy
    polynomial, x = 0, sp.symbols('x')
    for power, coeff in polynomial_dict.items():
        polynomial += x**power * coeff

    # construct remainder polynomial
    remainder_poly_dict = gen_remainder(values_list, polynomial_dict)


    return sp.Poly(polynomial, x)

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

def gen_remainder(values_list, sp_poly):
    coeff_dict = gen_poly_dict(sp_poly)
    remainder_dict = {}
    for key, coeff in coeff_dict.items():
        terms = coeff.args
        for term in terms:
           if term not in values_list:
               if key in remainder_dict:
                   remainder_dict[key] += term
              else:
                   remainder_dict[key] = term


    return remainder_dict


def factorize_polynomial(polynomial):
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
    
    # Print the factored polynomial
    sp.pprint(factored_poly, use_unicode=True)
    
    return factored_poly

# Generate the polynomial
polynomial = gen_polynomial(spin_S)
coeff_dict = print_polynomial(polynomial)
factored_poly = factorize_polynomial(polynomial)
