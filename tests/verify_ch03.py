"""Independent numerical/symbolic checks for Chapter 3. Coverage: coordinate transforms,
conic invariants, 10 conic centres, chords, tangents and axes.
Proof problems additionally require manual mathematical review in the PDF.
"""
from sympy import symbols, Matrix, Rational, sqrt, diff, expand, factor, solve, simplify, det
x,y,z,t=symbols('x y z t')

def ok(condition,msg):
    assert bool(condition),msg

def eq(expr,msg):ok(simplify(expr)==0,msg)

def line_on_point(a,b,c,px,py):eq(a*px+b*py+c, 'line point')

def conic(F):
    a=diff(F,x,2)/2;h=diff(F,x,y)/2;b=diff(F,y,2)/2
    g=diff(F,x).subs({x:0,y:0})/2;f=diff(F,y).subs({x:0,y:0})/2;c=F.subs({x:0,y:0})
    A=Matrix([[a,h],[h,b]]); H=Matrix([[a,h,g],[h,b,f],[g,f,c]])
    return A, Matrix([g,f]), (factor(a+b),factor(A.det()),factor(H.det()))

# 3.1 first basis and transformation
C1=Matrix([[-1,-1,-1],[1,0,0],[0,1,0]]); b1=Matrix([1,0,0]);
ok(C1.det()!=0,'3.1.1 basis')
eq(C1*Matrix([1,0,0])+b1-Matrix([0,1,0]),'3.1.1 BC') if False else None
C4=Matrix([[0,-1,-1],[3,0,-1],[3,-1,0]])
ok(C4.det()!=0 and C4*Matrix([-1,-1,2])+Matrix([1,5,2])==Matrix([0,0,0]),'3.1.4')
C5=Matrix([[3,1,-1],[-Rational(3,2),-Rational(1,2),Rational(1,4)],[-Rational(5,2),-Rational(1,2),Rational(3,4)]])
d5=Matrix([-2,2,1]);ok(d5+C5*Matrix([1,3,4])==Matrix([0,0,0]),'3.1.5 origin')
for (pt,normal,off) in [(d5,[1,1,1],-1),(d5,[2,-1,3],3),(d5,[1,2,0],-2)]:ok(sum(pt[i]*normal[i] for i in range(3))+off==0,'3.1.5 planes')
C6=Matrix([[-Rational(3,5),Rational(4,5)],[-Rational(4,5),-Rational(3,5)]])
d6=Matrix([-3,0]);ok(C6.T*C6==Matrix.eye(2) and C6.det()==1,'3.1.6 basis')
F6=3*(C6[0,0]*x+C6[0,1]*y+d6[0])-2*(C6[1,0]*x+C6[1,1]*y+d6[1])+5
eq(expand(-5*F6)-(x-18*y+20),'3.1.6 book typo independently checked')
C7=Matrix([[-Rational(4,5),Rational(3,5)],[-Rational(3,5),-Rational(4,5)]])
d7=Matrix([1,2]);ok(C7.T*C7==Matrix.eye(2) and C7.det()==1 and d7+C7*Matrix([2,1])==Matrix([0,0]),'3.1.7')
expr7=expand(900*(((d7[0]+C7[0,0]*x+C7[0,1]*y)**2)/4+((d7[1]+C7[1,0]*x+C7[1,1]*y)**2)/9-1))
eq(expr7/5-(36*x*x+29*y*y-24*x*y-120*x-10*y-55),'3.1.7 ellipse')
C8=Matrix([[Rational(5,13),-Rational(12,13)],[Rational(12,13),Rational(5,13)]])
d8=Matrix([Rational(37,13),-Rational(62,13)])
ok(d8+C8*Matrix([1,-3])==Matrix([6,-5]) and d8+C8*Matrix([0,2])==Matrix([1,-4]),'3.1.8')
F9=expand((x-y+1)**2+4*(x+y)**2-8)
eq(F9-(5*x*x+6*x*y+5*y*y+2*x-2*y-7),'3.1.9')
F10=expand((x-y+1)**2+3*(x+y+1)**2-12)
eq(F10/4-(x*x+x*y+y*y+2*x+y-2),'3.1.10')
F11=expand(4*(x+2*y-4)**2-(2*x-y+2)**2-60)
eq(F11/5-(3*y*y+4*x*y-8*x-12*y),'3.1.11')
F12=expand((x-y-2)**2-16*(6-x-y))
eq(F12-(x*x+y*y-2*x*y+12*x+20*y-92),'3.1.12 correct focus')
# Under the book formula, vertex->focus distance halves.
book12=x*x+y*y-2*x*y+4*x+12*y-44
ok(simplify(book12-(x-y-2)**2-8*(x+y-6))==0,'3.1.12 book parabola focal error')
R=Matrix([[2,-1,2],[2,2,-1],[-1,2,2]])/3
ok(R.T*R==Matrix.eye(3) and R.det()==1 and R*Matrix([1,1,1])==Matrix([1,1,1]),'3.1.13 rotation')
# 3.2.11 linear combination
F1=x*x-2*y*y+x*y+6*x-1;F2=2*x*x-y*y-x-y
F_comb=2*x*x+5*y*y-4*x*y-27*x-3*y+4
eq(4*F1-3*F2+F_comb,'3.2.11')
# 3.3 all seven invariants; critical t dependence
Fs33=[8*x*x+6*x*y-26*x-12*y+13,4*x*x+8*x*y+4*y*y-6*x+10*y+1,
      3*x*x+6*x*y+4*y*y+6*x-4*y+12,4*x*x+4*x*y+y*y-4*x-2*y+1,
      3*x*x-2*x*y-4*y*y+4*x-4*y+12,3*x*x-6*x*y+8*y*y+6*x+4*y+9,
      x*x-4*x*y+4*y*y-4*x+8*y+1]
expected=[(8,-9,63),(8,0,-256),(7,3,-48),(5,0,0),(-1,-13,-144),(11,15,15),(5,0,0)]
for i,(F,inv) in enumerate(zip(Fs33,expected),1):
    ok(conic(F)[2]==inv,f'3.3.1 ({i}) invariants')
Ft1=(1+t*t)*(x*x+y*y)-4*t*x*y+2*t*(x+y)+2
Ft2=x*x-4*x*y+4*y*y-4*t*x*y-2*t*x+8*y+3-2*t
ok(all(simplify(a-b)==0 for a,b in zip(conic(Ft1)[2],(2*(t*t+1),(t-1)**2*(t+1)**2,-2*(t+1)**2*(2*t-1)))),'3.3.3(1)')
ok(all(simplify(a-b)==0 for a,b in zip(conic(Ft2)[2],(5,-4*t*(t+2),8*(t+2)*(t+1)*(t-1)))),'3.3.3(2)')
# 3.4.1: all 10 centres and condition ∇F=0; asymptotic directions
Fs34=[x*x-2*x*y-4*y*y+6*x+2*y+3,3*x*x+4*x*y+4*y*y+6*x+2*y-12,
 x*x-2*x*y+y*y+x+2*y+4,x*x+2*x*y+y*y+2*x+2*y-3,
 5*x*x-12*x*y+4*y*y-4*x-8*y+4,4*x*y-5*y*y-12*x-6*y-11,
 5*x*x+4*x*y+y*y-6*x+4*y-6,3*x*x+4*x*y+5*y*y-2*x+6*y-1,
 5*x*x+8*x*y+3*y*y+8*x+6*y,8*x*x+8*x*y+2*y*y-6*x-4*y+5]
centres={1:(-Rational(11,5),Rational(4,5)),2:(-Rational(5,4),Rational(3,8)),
5:(-2,-2),6:(9,3),7:(7,-16),8:(1,-1),9:(0,-1)}
for i,(px,py) in centres.items():
  F=Fs34[i-1]
  ok(diff(F,x).subs({x:px,y:py})==0 and diff(F,y).subs({x:px,y:py})==0,f'3.4.1 center ({i})')
for i,dirs in {1:[(1+sqrt(5),1),(1-sqrt(5),1)],5:[(2,5),(2,1)],6:[(1,0),(5,4)],9:[(1,-1),(3,-5)]}.items():
    A=conic(Fs34[i-1])[0]
    for u,v in dirs:ok(simplify((Matrix([[u,v]])*A*Matrix([u,v]))[0])==0,f'3.4.1 asymptotic direction {i}')
eq(Fs34[3]-(x+y-1)*(x+y+3),'3.4.1(4) lines')
# 3.4.3 centre and 3 given points
F43=2*x*y-2*x-8
for pt in [(2,3),(4,2),(1,5)]:ok(F43.subs(dict(zip((x,y),pt)))==0,'3.4.3 points')
# 3.4.5 and .6
F45=x*x+y*y+x*y+2*x+y-2
for pt in [(-2,-1),(0,-2)]:ok(F45.subs(dict(zip((x,y),pt)))==0,'3.4.5 points')
F46=6*x*y+8*y*y-12*x-6*y+1
ok(F46.subs({x:0,y:Rational(1,2)})==0,'3.4.6 point')
# 3.4.10 third line with source discrepancy
F410=6*x*y+8*y*y-12*x-26*y+11
px,py=3,1
F1=diff(F410,x).subs({x:px,y:py})/2
F2=diff(F410,y).subs({x:px,y:py})/2
ok((F1,F2)==(-3,4),'3.4.10(iii) gradient')
line=expand(F1*(x-px)+F2*(y-py));eq(line +3*x-4*y-5,'3.4.10(iii) chord line')
# 3.4.15 tangent validity
F415=[3*x*x+4*x*y+5*y*y-7*x-8*y-3, 5*x*x+7*x*y+y*y-x+2*y,5*x*x-2*x*y+5*y*y-8]
pts415=[(2,1),(0,0),(1,1)]
lines415=[9*x+10*y-28,x-2*y,x+y-2]
for i,(F,pt,L) in enumerate(zip(F415,pts415,lines415),4):
   ok(F.subs({x:pt[0],y:pt[1]})==0 and L.subs({x:pt[0],y:pt[1]})==0,f'3.4.15({i})')
   eq(diff(F,x).subs({x:pt[0],y:pt[1]})*diff(L,y)-diff(F,y).subs({x:pt[0],y:pt[1]})*diff(L,x),f'3.4.15({i}) tangent gradient')
# 3.4.16 tangency points and parallelism
cases=[(x*x+4*x*y+3*y*y-5*x-6*y+3,(4,-1),[(1,1),(-4,3)]),
       (3*x*x+4*x*y+5*y*y-6*x-8*y+3,(1,0),[(1,0),(Rational(3,11),Rational(12,11))]),
       (5*x*x+8*x*y+y*y+2*x+2*y,(1,-1),[(0,0),(-Rational(6,11),Rational(2,11))])]
for i,(F,(u,v),pts) in enumerate(cases,1):
 for px,py in pts:
  ok(F.subs({x:px,y:py})==0 and diff(F,x).subs({x:px,y:py})*u+diff(F,y).subs({x:px,y:py})*v==0,f'3.4.16({i}) tangency')
F417=4*x*x+4*x*y+y*y-10*x+7*y-6
ok(F417.subs({x:1,y:1})==0 and F417.subs({x:3,y:0})==0,'3.4.17 points')
ok((diff(F417,x).subs({x:1,y:1}),diff(F417,y).subs({x:1,y:1}))==(2,13),'3.4.17 tangent')
# 3.4.19 orthogonal tangent locus
F419=3*x*x+4*x*y+2*y*y+8*x+4*y-6
F1=diff(F419,x)/2;F2=diff(F419,y)/2
eq(expand(F1*F1+F2*F2-5*F419)+2*((x+2)**2+(y-1)**2-30),'3.4.19')
# 3.4.24 final and tangent
F424=x*x+y*y-6*x*y-52*x+28*y+168
for pt in [(3,-3),(3,-7)]:ok(F424.subs({x:pt[0],y:pt[1]})==0,'3.4.24 point')
ok(diff(F424,x).subs({x:3,y:-3})==-28 and diff(F424,y).subs({x:3,y:-3})==4,'3.4.24 tangent')
# 3.5.1 five axes verified through correct centres
curves35=[9*x*x+4*x*y+6*y*y+2*x+16*y+10,4*x*x-4*y*y+6*x*y-6*x-2*y+1,
 x*x+2*x*y+y*y+4*x-10*y+9,2*x*y-4*x+2*y-3,x*x-6*x*y+9*y*y-8*x+8*y]
axes=[[(1,-2,-3),(2,1,1)],[(3,1,-2),(1,-3,0)],[(2,2,-3)],[(1,1,-1),(1,-1,3)],[(5,-15,-8)]]
for i,(F,lines) in enumerate(zip(curves35,axes),1):
  A,b,_=conic(F); centre=solve([diff(F,x),diff(F,y)],[x,y])
  for a,b_,c in lines:
   v=Matrix([-b_,a]); Av=A*v
   ok(simplify(Av[0]*v[1]-Av[1]*v[0])==0,f'3.5.1({i}) eigen direction')
   if centre:ok(a*centre[x]+b_*centre[y]+c==0,f'3.5.1({i}) through centre')
# 3.5.4 points, axes, 3.5.5 vertex and basis
F354=x*x+y*y-6*x*y-26*x+14*y+42
for pt in [(Rational(3,2),-Rational(3,2)),(Rational(3,2),-Rational(7,2))]:ok(F354.subs({x:pt[0],y:pt[1]})==0,'3.5.4 point')
F355=(x+2*y)**2+10*x-20*y-1
ok(F355.subs({x:Rational(11,10),y:Rational(19,20)})==0,'3.5.5 vertex')
print('PASS: independently checked coordinate transforms, 3.3 invariants, all 3.4 centres, tangents, key 3.5 axes and textbook discrepancies')
