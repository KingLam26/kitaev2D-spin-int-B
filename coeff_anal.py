##### spin 1.5 #####
print("\n>>>spin 1.5\n")

d1 = 22973.4985557602266619427492945
d2 = 28593.7042944068864775129139834
d3 = 35920.9438067482684683850866605
d4 = 35970.0067626467559813040863878
d5 = 33459.8776620620398467128998918
d6 = 42250.1284261502151881160627356
d7 = 49770.2872895085381091211911294

print((2*d2+2*d6) / d3)
print((d1+d7+2*d5+2*d4) / d3)

##### spin 2 #####
print("\n>>>spin 2\n")

d1 = 1.8604521
d2 = 4.1379024
d3 = 9.2064512

d3_d2 = d3/d2
d2_d1 = d2/d1

e11 = d3/d2 - d2_d1

print(f"d2_d1: {d2_d1}, d3_d2: {d3_d2}")
print(f"e11: {e11}")

##### spin 3 #####
print("\n>>>spin 3\n")

d1 = 36052814083126422740
d2 = 124066871221800233745
d3 = 427101825544440584157

d3_d2 = d3/d2
d2_d1 = d2/d1

e11 = d3/d2 - d2_d1

print(f"d2_d1: {d2_d1}, d3_d2: {d3_d2}")
print(f"e11: {e11}")

##### spin 4 #####
print("\n>>>spin 4\n")

d1 = 1.290249294643517
d2 = 6.0074698994905354
d3 = 9.480000447372665
d4 = 27.98033970269736
d5 = 44.15883619644383
d6 = 69.69445768607373

a = d2/d1
b = d3/d1
print(f"a: {a} b: {b}")

e_aa = (d4/d1) - a**2
e_ab = (d5/d1) - a*b
e_bb = (d6/d1) - b**2
print(f"e_aa: {e_aa} e_ab: {e_ab} e_bb: {e_bb}")

e_a1b1 = e_ab / e_aa
e_b1b1 = e_bb / e_aa
print(f"e_a1b1: {e_a1b1} e_b1b1: {e_b1b1}")

delta = (e_b1b1/e_a1b1) - e_a1b1
print(f"delta: {delta}")


from sympy import symbols

# Define variables with subscripts
a = symbols('a0:3')
x = symbols('x')

poly = x*a[0]

# Print the variables
print(poly)
