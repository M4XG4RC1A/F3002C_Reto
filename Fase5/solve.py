from sympy import sqrt
from sympy import symbols, Eq, solve

# Se definen las variables
wr, ws, wi, V_plus, Va, Vb = symbols('wr ws wi V_plus Va Vb')
w0 = symbols('w0')

# A quitarle el trabajo a mathematica
eq1_new = Eq(V_plus, 1/sqrt(3) * (wr + ws + wi - 3*w0))
eq2_new = Eq(Va, 1/2 * (1 - 1/sqrt(3))*wr + 1/2 * (-1 - 1/sqrt(3))*ws + 1/sqrt(3)*wi)
eq3_new = Eq(Vb, 1/2 * (1 + 1/sqrt(3))*wr + 1/2 * (-1 + 1/sqrt(3))*ws - 1/sqrt(3)*wi)
solutions_new = solve((eq1_new, eq2_new, eq3_new), (wr, ws, wi))
print(solutions_new)
