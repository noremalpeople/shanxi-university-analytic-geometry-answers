#!/usr/bin/env python3
"""Independent SymPy spot checks for Exercise 1.1.

These tests verify the algebraic/computational examples, not all prose proofs.
"""
from sympy import Matrix, Rational, symbols, simplify, solve, det, factor


def eq(a, b):
    """Identity of scalar or matrix under simplification."""
    d = a - b
    if isinstance(d, Matrix):
        return all(simplify(x) == 0 for x in d)
    return simplify(d) == 0


def vec(x, y):
    return Matrix([x, y])


def cross2(a, b):
    return factor(a[0] * b[1] - a[1] * b[0])


def ratio(a, b, p):
    """Signed ratio for collinear points in 2D."""
    d = b - a
    k = 0 if d[0] != 0 else 1
    t = factor((p[k] - a[k]) / d[k])
    assert eq(p - a, t * d), f"Not collinear: {a}, {b}, {p}"
    return factor(t / (1 - t))


def test():
    alpha, beta, x, y = symbols("alpha beta x y", real=True)
    # 1: α = AB + BC, β = BC - AB
    s1 = solve((x + y - alpha, -x + y - beta), (x, y))
    assert s1[x] == (alpha - beta) / 2
    assert s1[y] == (alpha + beta) / 2
    print("[PASS] 01 parallelogram diagonals")

    # 2: AE=AB+AD/2 and AF=AB/2+AD
    s2 = solve((x + y / 2 - alpha, x / 2 + y - beta), (x, y))
    assert eq(s2[x], (4 * alpha - 2 * beta) / 3)
    assert eq(s2[y], (4 * beta - 2 * alpha) / 3)
    print("[PASS] 02 parallelogram midpoints")

    # 3, 4: basis vector checks
    a, d = vec(2, 0), vec(0, 3)
    C, B = d + a / 2, a
    E, F = (B + C) / 2, (C + d) / 2
    assert eq(E, 3 * a / 4 + d / 2)
    assert eq(F, a / 4 + d)
    assert eq(a, 4 * (2 * E - F) / 5)
    assert eq(C, 2 * (E + 2 * F) / 5)
    assert eq(C-B, 2 * (4 * F - 3 * E) / 5)
    assert eq(d-B, 2 * (F-E))
    print("[PASS] 03-04 triangle and trapezoid vector relations")

    # 5: direct contradiction between printed input and end-of-book answer
    a, b = vec(1, 0), vec(0, 1)
    A, B, C = vec(0, 0), a, a + b
    D_key, E_key, F_key = 2 * b, -a + 2 * b, b - a
    assert eq(E_key - D_key, -a)
    assert not eq(E_key - D_key, a)
    assert eq(C-F_key, 2*a)
    assert cross2(C - D_key, F_key - A) == 0
    assert cross2(B - C, E_key - F_key) == 0
    # Formally retaining DE=a forces a different answer.
    D_orig, E_orig = -2*a + 4*b, -a + 4*b
    assert eq(E_orig-D_orig, a)
    assert eq(C-F_key, 2*a)
    assert eq(D_orig-C, 3*(b-a))
    assert eq(E_orig-C, 3*b-2*a)
    assert cross2(D_orig-C, A-F_key) == 0
    assert cross2(E_orig-F_key, C-B) == 0
    print("[PASS] 05 confirmed end-answer conflicts with scanned condition")

    # 6: arbitrary midpoint identity
    A,B,C,D=vec(1,4),vec(5,-1),vec(-2,3),vec(1,9)
    P,Q=(A+B)/2,(C+D)/2
    assert eq(2*(Q-P),(C-A)+(D-B))
    print("[PASS] 06 midpoint identity")

    # 7-14: centroid grouping and assignment invariance
    pts=[vec(1,0),vec(3,2),vec(6,-1),vec(-2,4),vec(9,-7)]
    M=sum(pts[:3],Matrix([0,0]))/3
    N=sum(pts[3:],Matrix([0,0]))/2
    G=sum(pts,Matrix([0,0]))/5
    assert eq(G,(3*M+2*N)/5)
    assert eq(G,M+Rational(2,5)*(N-M))
    dest=[vec(1,1),vec(2,-4),vec(-5,-8),vec(10,6),vec(3,3)]
    lhs=sum((dest[i]-pts[i] for i in range(5)),Matrix([0,0]))
    permutation=[3,0,4,2,1]
    rhs=sum((dest[permutation[i]]-pts[i] for i in range(5)),Matrix([0,0]))
    assert eq(lhs,rhs)
    print("[PASS] 07-14 centroid and matching sum")

    # 18: directed ratio identities
    A,B,C=vec(0,0),vec(1,0),vec(Rational(5,2),0)
    assert eq(ratio(A,B,C)*ratio(B,A,C),1)
    assert eq(ratio(A,B,C)+ratio(A,C,B),-1)
    print("[PASS] 18 oriented ratios")

    # 19: independent 2D line intersection
    A,B,C=vec(0,0),vec(1,0),vec(0,1)
    F=vec(Rational(2,3),0)
    E=vec(0,Rational(1,3))
    D=vec(Rational(3,4),Rational(1,4))
    t,u=symbols('t u')
    s19=solve(list(A+t*(D-A)-(E+u*(F-E))),(t,u),dict=True)[0]
    G=A+s19[t]*(D-A)
    assert ratio(E,F,G)==Rational(3,2)
    assert ratio(A,D,G)==Rational(8,7)
    print("[PASS] 19 ratios 3/2, 8/7")

    # 20: four listed ratios
    F=vec(Rational(1,4),0)
    O=C+Rational(2,3)*(F-C)
    D=vec(Rational(1,3),Rational(2,3))
    E=vec(0,Rational(2,5))
    assert [ratio(A,D,O),ratio(B,C,D),ratio(C,A,E),ratio(B,E,O)] == [1,2,Rational(3,2),5]
    print("[PASS] 20 ratios 1, 2, 3/2, 5")

    # 21: independently verify both barycentric-directed-ratio formulae.
    mu,nu=symbols('mu nu')
    lam=1-mu-nu
    A,B,C=vec(0,0),vec(1,0),vec(0,1)
    E=(mu*B+nu*C)/(mu+nu)
    D=mu*B+nu*C
    assert eq(D,(1-lam)*E)
    assert eq((1-lam)/lam,(mu+nu)/lam)
    assert eq(nu/(mu+nu)/(1-nu/(mu+nu)),nu/mu)
    print("[PASS] 21 barycentric coefficients and ratios")

    # 22: product 1 without Euclidean concurrency due parallel lines
    A,B,C=vec(0,0),vec(1,0),vec(0,1)
    D,E,F=vec(2,0),vec(2,-1),vec(0,Rational(1,2))
    ratios=[ratio(A,B,D),ratio(B,C,E),ratio(C,A,F)]
    assert ratios==[-2,-Rational(1,2),1]
    assert factor(ratios[0]*ratios[1]*ratios[2])==1
    assert cross2(E-A,F-B)==0 and cross2(E-A,D-C)==0
    assert cross2(E-A,B-A)!=0  # AE and BF do not coincide
    print("[PASS] 22 parallel Ceva counterexample")

    # 23: symbolic 3-point collinearity polynomial equals negative-product equation.
    p,q,r=symbols('p q r')
    lhs=det(Matrix([[p,0,1],[1-q,q,1],[0,1-r,1]]))
    assert eq(lhs,p*q*r+(1-p)*(1-q)*(1-r))
    print("[PASS] 23 Menelaus determinant and signed-product identity")


if __name__ == "__main__":
    test()
    print("All computable spot checks PASSED. This is not a machine proof of every prose step.")
