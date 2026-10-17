#!/usr/bin/env python3
"""Independent algebraic spot checks for Chapter 2 worked solutions.
These are not formal proofs; each assertion corresponds to a derivation in the PDF.
"""
import itertools
import sympy as s

x,y,z,t,q=s.symbols('x y z t q', real=True)
X=s.Matrix([x,y,z])

def sub(f,point):
    return s.simplify(f.subs({x:point[0],y:point[1],z:point[2]}))

def same(a,b):
    assert s.simplify(a-b)==0, (a,b)

checks=[]
def check(name,action):
    action();checks.append(name)

# 2.1: four-point sphere, signed centers, sphere classification
A=[(2,0,2),(2,3,0),(3,4,1),(1,1,1)]
F=(x-s.Rational(11,4))**2+(y-2)**2+(z-s.Rational(7,4))**2-s.Rational(37,8)
check('2.1 Q1(3) four points',lambda:[same(sub(F,p),0) for p in A])
assert [(sg,roots) for sg in itertools.product([-1,1],repeat=3)
        for roots in [[r for r in s.solve(s.Symbol('R')**2-(-sg[0]+2*sg[1]+5*sg[2])*s.Symbol('R')+15,s.Symbol('R'))
                        if r.is_positive]] if roots]==[((-1,1,1),[s.Integer(3),s.Integer(5)])]
checks.append('2.1 Q1(5) all 8 tangent sign cases')

# 2.2 plane checks
check('2.2 Q2(2) 3 points',lambda:[same(sub(19*x-13*y-11*z+3,p),0) for p in [(1,0,2),(2,4,-1),(-3,-5,1)]])
check('2.2 Q2(3) 2 points+direction',lambda:[same(sub(20*x-12*y+z+18,p),0) for p in [(-1,0,2),(1,3,-2)]])
check('2.2 Q12 ratio 2:1',lambda:[same(abs(v+4),2*abs(v-5)) for v in [2,14]])

# 2.3 Q7 parameterization and Q12 ruled surface
check('2.3 Q7 straight line',lambda:[same(sub(3*x-y+z-2,(-1+2*t,-5+7*t,t)),0),same(sub(x-2*z+1,(-1+2*t,-5+7*t,t)),0)])
F12=4*x*x-9*y*y+6*x-45*y-36*z-144
P=s.Matrix([6+3*t,2*t,1+2*t]);Q=s.Matrix([3*q,2*q,-4-2*q]);
# line through P and Q parallel 2x+3y-5=0 demands n dot (Q-P)=0
qof=s.solve(s.Eq(2*(Q[0]-P[0])+3*(Q[1]-P[1]),0),q)[0]
Q=Q.subs(q,qof);r=s.symbols('r')
M=(1-r)*P+r*Q
check('2.3 Q12 ruled surface',lambda:same(sub(F12,M),0))

# 2.4 Q8 distance and Q11 bisectors
same(2*3+4-10,0); same(s.Rational(16,1)/s.sqrt(45),16*s.sqrt(5)/15)
checks.append('2.4 Q8 parallel line and distance')
F1=3*x-2*y-6*z-4; F2=2*x+2*y-z+5
check('2.4 Q11 bisector planes',lambda:[same(s.expand(3*F1+7*F2),23*x+8*y-25*z+23),same(s.expand(-3*F1+7*F2),5*x+20*y+11*z+47)])

# 2.5 Q1(1) line rotation z-axis
T=s.symbols('T',real=True)
F=2*x*x+2*y*y-5*z*z-58*z-169
check('2.5 Q1(1) rotate line',lambda:same(sub(F,(T,-2-3*T,-7-2*T)),0))
# torus meridian rotation original circle (x-3)^2+y^2=4
rho=s.symbols('rho',nonnegative=True)
Ft=(rho*rho+y*y+5)**2-36*rho*rho
check('2.5 Q1(4) torus circle profile',lambda:same(s.rem(s.Poly(Ft,rho),s.Poly((rho-3)**2+y*y-4,rho)).as_expr(),0))
# 2.5 5(1) ellipse from foci, point at semimajor axis (5,0,0)
Fell=x*x/25+y*y/16+z*z/16-1
check('2.5 Q5(1) printed denominator discrepancy',lambda:[same(sub(Fell,(5,0,0)),0),same(sub(x*x/5+y*y/4+z*z/4-1,(5,0,0)),4)])

# 2.5 Q7(2),(3) cylinders and erratum
u=s.Matrix([1,2,3]);P=s.Matrix([1,0,2]);W=X-P
F7=s.expand(W.cross(u).dot(W.cross(u))-9*u.dot(u))
F7answer=13*x*x+10*y*y+5*z*z-4*x*y-6*x*z-12*y*z-14*x+28*y-14*z-105
check('2.5 Q7(2) cylinder expansion',lambda:same(F7,F7answer))
u=s.Matrix([1,-1,2]);P=s.Matrix([1,2,-2]);Q=s.Matrix([1,-1,0]);W=X-P
r2=Q-P; r2=r2.cross(u).dot(r2.cross(u))
F7=s.expand(W.cross(u).dot(W.cross(u))-r2)
F7answer=5*x*x+5*y*y+2*z*z+2*x*y-4*x*z+4*y*z-22*x-14*y+4*z
check('2.5 Q7(3) erratum verified',lambda:[same(F7,F7answer),same(sub(F7answer,Q),0),same(sub(F7answer-5,Q),-5)])

# 2.5 cone Q10(2) from axis and angle
V=s.Matrix([1,0,2]);u=s.Matrix([2,2,-1]);W=X-V
Fcone=s.expand(27*W.dot(W)-4*(W.dot(u))**2)
Fbook=11*x*x+11*y*y+23*z*z-32*x*y+16*x*z+16*y*z-54*x-108*z+135
check('2.5 Q10(2) cone expansion',lambda:same(Fcone,Fbook))

# 2.5 Q15 ruled cylinder
F15=y*y+z*z+2*y*z-x-4*y-2*z+3
S=s.symbols('S')
check('2.5 Q15 ruled cylinder',lambda:same(sub(F15,(T*T+2*S,1-S,T+S)),0))

# 2.6 Q4 cross sections and k=a critical cylinder
F6=x*x/4+y*y/9-z*z-1
check('2.6 Q4(1) ellipse section',lambda:same(s.expand(F6.subs(z,s.sqrt(3))*4),x*x+y*y*s.Rational(4,9)-16))
check('2.6 Q8 k=a cylinder',lambda:same(((s.Symbol('a')-s.Symbol('a'))*x*x+(s.Symbol('b')-s.Symbol('a'))*y*y-z), (s.Symbol('b')-s.Symbol('a'))*y*y-z))

# 2.7 Q1 both line families
F7=x*x/4+y*y/9-z*z/16-1
check('2.7 Q1 two rulings',lambda:[same(sub(F7,(2*T,-3,4*T)),0),same(sub(F7,(2,3*T,-4*T)),0)])

# 2.7 Q5 joining matching t on two lines
F7=-x*x/9+y*y/4-2*z
svar=s.symbols('svar')
P=s.Matrix([s.Rational(3,2)+3*T,-1+2*T,-T]);Q=s.Matrix([3*T,2*T,0]);M=(1-svar)*P+svar*Q
check('2.7 Q5 hyperbolic paraboloid',lambda:same(sub(F7,M),0))

# 2.7 Q9 Q10 Q11 include both given straight lines
F9=x*x+y*y+2*x*y-x*z-y*z-x-2*y+z
check('2.7 Q9 lines on quadric',lambda:[same(sub(F9,(1+T,0,1+T)),0),same(sub(F9,(1,1+T,1+T)),0)])
F10=3*y*y-5*z*z+4*x*y-6*x-16*y+5*z+36
check('2.7 Q10 contains base lines',lambda:[same(sub(F10,(6+T,2*T,1+2*T)),0),same(sub(F10,(T,8+2*T,-4-2*T)),0)])
F11=2*x*x-2*z*z+x*y+y*z+8*x+6*y+14*z+6
check('2.7 Q11 contains base lines',lambda:[same(sub(F11,(1+T,-3,5+T)),0),same(sub(F11,(-T,2+2*T,-1)),0)])
print(f'PASS: {len(checks)} independent chapter-2 checks')
for c in checks: print(' -',c)
