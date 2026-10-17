#!/usr/bin/env python3
"""Independent exact symbolic/numeric regression checks for exercises 1.2-1.5."""
from sympy import Matrix as M, Rational as Q, symbols, simplify, sqrt, cos, pi, det

def zero(v):
    if isinstance(v, M):
        return all(simplify(x)==0 for x in v)
    return simplify(v)==0

def T(a,b,c):
    return a.dot(b.cross(c))

def check():
    # 1.2 affine-basis transform, part 1: parallelogram diagonals
    B,D,A,C=M([1,0]),M([0,1]),M([0,0]),M([1,1])
    P=(B+4*D)/5; R=(A+5*C)/6
    assert P==M([Q(1,5),Q(4,5)]) and R-P==M([Q(19,30),Q(1,30)])
    e=M([1,Q(1,2)]);f=M([Q(1,2),1]);N=M.hstack(e,f)
    assert N.inv()*M([1,0])==M([Q(4,3),-Q(2,3)])
    assert N.inv()*M([0,1])==M([-Q(2,3),Q(4,3)])
    N=M.hstack(M([Q(2,3),Q(1,2)]),M([Q(1,6),1]))
    for old,new in [(M([-Q(1,3),-1]),M([-Q(2,7),-Q(6,7)])),(M([Q(2,3),-1]),M([Q(10,7),-Q(12,7)])),(M([-Q(1,3),0]),M([-Q(4,7),Q(2,7)]))]:
        assert N.inv()*old==new
    V=M.hstack(M([2,1]),M([1,2]))
    assert V.inv()*M([1,0])==M([Q(2,3),-Q(1,3)])
    assert V.inv()*M([0,1])==M([-Q(1,3),Q(2,3)])
    a,b,c=M([3,7,3]),M([10,0,3]),M([8,2,3]); assert zero(7*c-2*a-5*b)
    A,B,C=M([3,4,1]),M([2,5,0]),M([6,1,4]);assert zero(C-A+3*(B-A))
    print('[PASS] 1.2: Q1-6, Q10')

    # 1.3 dot-product concrete questions
    a2,b2,ab=9,4,3
    assert 9*a2+12*ab+4*b2==133
    assert 4*a2-20*ab+25*b2==76
    assert 6*a2-11*ab-10*b2==-19
    assert simplify((25+12*sqrt(2))/9-(9+16+12*sqrt(2))/9)==0
    a,b=M([-1,2]),M([1,1]);pr=(a.dot(b)/b.dot(b))*b;assert pr==M([Q(1,2),Q(1,2)])
    a,b,c=M([1,0,2]),M([1,1,1]),M([1,0,-2]);n=b.cross(c)
    p2=a.dot(n)/n.dot(n)*n; p1=a-p2
    assert p2==M([Q(4,7),-Q(6,7),Q(2,7)]) and p1==M([Q(3,7),Q(6,7),Q(12,7)])
    x,y,z=M([3,-5,7]),M([1,2,-9]),M([-4,2,2]);d=-(x+y+z)
    assert zero(x.dot(x)+y.dot(y)-z.dot(z)-d.dot(d)+2*(x+z).dot(y+z))
    assert zero(x.dot(x)+z.dot(z)-y.dot(y)-d.dot(d)+2*(x+y).dot(y+z))
    assert zero(y.dot(y)+z.dot(z)-x.dot(x)-d.dot(d)+2*(x+y).dot(x+z))
    print('[PASS] 1.3: Q1, Q3-5, Q12')

    # 1.4 and 1.5 known determinant values
    aa,bb,cc=M([1,1,1]),M([6,-1,2]),M([5,1,7]); assert (bb-aa).cross(cc-aa).dot((bb-aa).cross(cc-aa))==884
    a,b,c=M([1,0,-1]),M([1,-2,0]),M([-1,2,1]);assert (3*a+b-c).cross(a-b+c)==M([16,4,16])
    a,b,c=M([3,2,-1]),M([0,3,3]),M([5,2,1]);assert T(a,b,c)==36
    assert -T(a,b,c)==-36 # left-handed coordinate physical scalar triple
    a,b,c=M([1,0,3]),M([2,4,3]),M([3,4,1]);assert T(a,b,c)==-20
    p=M([1,0,1]);a=M([-1,1,5])-p;b=M([-1,-3,-3])-p;c=M([0,3,4])-p
    assert abs(T(a,b,c))==32
    triples=[(M([1,2,3]),M([2,3,4]),M([3,4,5]),0),
             (M([3,2,-1]),M([0,3,7]),M([2,2,1]),1),
             (M([2,1,3]),M([1,-1,-4]),M([5,1,2]),0)]
    for a,b,c,expected in triples:assert T(a,b,c)==expected
    print('[PASS] 1.4 Q2-3, 1.5 Q4-6')

    # 1.5 Q2 exact derivation and printed-denominator counterexample
    a,b,c=M([2,0,0]),M([0,2,0]),2
    xi=(c*a-a.cross(b))/a.dot(a)
    assert xi==M([1,0,-1]) and a.dot(xi)==c and a.cross(xi)==b
    xifalse=(c*a-a.cross(b))/sqrt(a.dot(a))
    assert a.dot(xifalse)!=c
    # 1.5 Q3 a dot beta nonzero, gamma perpendicular beta
    alpha,beta,gamma=M([1,1,0]),M([0,2,0]),M([3,0,1]);val=5
    xi0=gamma.cross(beta)/beta.dot(beta)
    s=(val-alpha.dot(xi0))/alpha.dot(beta)
    xi=xi0+s*beta
    assert alpha.dot(xi)==val and beta.cross(xi)==gamma
    print('[PASS] 1.5 Q2 corrected squared denominator, Q3')

    # Verify all named vector identities for several arbitrary vectors.
    data=[(M([2,-3,1]),M([4,1,-1]),M([-2,5,7]),M([6,-4,3])),
          (M([0,0,0]),M([2,2,-1]),M([0,1,3]),M([-2,-2,4])),
          (M([1,0,0]),M([0,1,0]),M([0,0,1]),M([1,2,3]))]
    for a,b,c,d in data:
        assert zero(a.cross(b.cross(c))+b.cross(c.cross(a))+c.cross(a.cross(b)))
        assert zero(a.cross(b.cross(c.cross(d)))-b.dot(d)*a.cross(c)+b.dot(c)*a.cross(d))
        assert zero(a.cross(b.cross(c.cross(d)))-T(a,c,d)*b+a.dot(b)*c.cross(d))
        assert zero(T(b,c,d)*a-a.dot(b)*c.cross(d)-a.dot(c)*d.cross(b)-a.dot(d)*b.cross(c))
        assert zero(a.cross(b).cross(c.cross(d))-T(a,b,d)*c+T(a,b,c)*d)
        assert zero(a.cross(b).cross(c.cross(d))-T(a,c,d)*b+T(b,c,d)*a)
        assert zero(T(b,c,d)*a-T(a,c,d)*b+T(a,b,d)*c-T(a,b,c)*d)
        assert zero(T(a.cross(b),b.cross(c),c.cross(a))-T(a,b,c)**2)
        assert zero(a.cross(b)+b.cross(c)+c.cross(a)-(a-c).cross(b-c))
        assert zero((a.cross(b)).dot(c.cross(d))+(a.cross(d)).dot(b.cross(c))+(a.cross(c)).dot(d.cross(b)))
    print('[PASS] 1.5 Q1 and Q7-14 evaluated on independent vector cases')

if __name__=='__main__':check()
