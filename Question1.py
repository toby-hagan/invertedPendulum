import sympy as sym

# step 0. introduce the necessary symbols
M, m, g, x1, x2, x3, x4, F, ell = sym.symbols('M, m, g, x1, x2, x3, x4, F, ell')
phi = 4*m*ell*x4**2*sym.sin(x3) + 4*F - 3*m*g*sym.sin(x3)*sym.cos(x3)
phi /= 4*(M+m) - 3*m*sym.cos(x3)**2

# declaration of psi function
psi = -3*(m*ell*x4**2*sym.sin(x3)*sym.cos(x3) + F*sym.cos(x3) - (M + m)*g*sym.sin(x3))
psi /= (4*(M + m) - 3*m*sym.cos(x3)**2)*ell

# step 1. differentiate phi wrt F, x3 and x4
d_phi_F = phi.diff(F)
d_phi_x3 = phi.diff(x3)
d_phi_x4 = phi.diff(x4)

# psi differentiations
d_psi_F = psi.diff(F)
d_psi_x3 = psi.diff(x3)
d_psi_x4 = psi.diff(x4)


# (question 1 )
# step 2. substitute F=0, x3=0, x4=0 into the
#         above derivatives
d_phi_F_eq = d_phi_F.subs([(F, 0), (x3, 0), (x4, 0)])
d_phi_x3_eq = d_phi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
d_phi_x4_eq = d_phi_x4.subs([(F, 0), (x3, 0), (x4, 0)])

d_psi_F_eq = d_psi_F.subs([(F, 0), (x3, 0), (x4, 0)])
d_psi_x3_eq = d_psi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
d_psi_x4_eq = d_psi_x4.subs([(F, 0), (x3, 0), (x4, 0)])

sym.pprint(d_phi_F_eq)
'''
sym.pprint(d_phi_x3_eq)
sym.pprint(d_phi_x4_eq)

sym.pprint(d_psi_F_eq)
sym.pprint(d_psi_x3_eq)
sym.pprint(d_psi_x4_eq)
'''
