# -*- coding: utf-8 -*-
"""
    The purpose of this code is to calculate supported forces on the axles of
    a tandem axle 40ft trailor. The calculation assumes that the trailor is a 
    beam of length l supported at both ends. First ends will be a fixed support (A). 
    Seconnd end will be a roller support (C). There is additional roller support offset
    from the end roller support (B) by the dimension a toward the center of the beam

    @author: yjeong
"""

from sympy import *

x = symbols ('x') #initialize X this is needed for integration
#w = weight per length (1000lbs/ft)
#l = length of the beam/trailor (40ft)
#a = space between tandem axles (1ft)
#c = constants generated during integration
w, l, a, c, Ma = symbols ('w l a c Ma') #initialize geometric constants
A, B, C = symbols ('A B C') #initialize reactionary forces and moments

#print('Solving Statically Indeterminate Beam using Double Integration Piece Wise Fx','\n')
#-----------------------------FIRST SECTION------------------------------------
#set Moment A equation from Sum of Moment at C
ma = A*l + B*a - w*l**2/2
mexp1 = ma + A*x-w*x**2/2

#print('From sum of static some of moment: ', ma)
#print('From derivitive of first piece wise: EIv"(x) = ',mexp1,'\n')
#print('Solving First Section',)
#----------------------------------FIRST INTEGRATION---------------------------
mexpint1 = integrate(mexp1,x)
mexpint1 = mexpint1+mexpint1.subs(x,0) #+c, which is 0 in this case
#print("EIv'(x) =  ", mexpint1)

#---------------------------------SECOND INTEGRATION--------------------------
mexpint2 = integrate(mexpint1,x)
mexpint2 = mexpint2+mexpint2.subs(x, 0) #+c, which is 0 in this case
#if c in the second integration is equal to itself
#print('EIv(x) = ', mexpint2,'\n')

#---------------SOLVING FOR A in respect to B and w ---------------------------
#print('Solve for A at L-a where deflection is 0') #above first axle
#print('EIv(l-a) = 0 = ',mexpint2.subs(x,(l-a)), '\n')
#print('sanity: ',expand(mexpint2.subs(x,(l-a))),'\n')
mexpint2 = expand(mexpint2)
eqA = simplify(mexpint2.subs([(A,0),(x,(l-a))])/(-1*mexpint2.subs([(x,(l-a)),(B,0),(w,0),(A,1)])))
#print('A = ',eqA,'\n')

#-------------------------------THIRD SECTION----------------------------------
mexp2 = w*x**2/2 - C*x #initialize equation
#print('Solving Third Section')
#print('Eiv"(x) = ',mexp2)
#----------------------------------FIRST INTEGRATION---------------------------
mexp2int1 = integrate(mexp2,x)
mexp2int1 = mexp2int1 + mexp2int1.subs(x,l)
#print("EIv'(x) =  ", mexpint1)
#mexp2c3 = -1*(simplify(mexp2int1.subs(x,(l)))) #since c3 != 0
#print(mexp2int1+mexp2c3)
#----------------------------------SECOND INTEGRATION--------------------------
mexp2int2 = integrate(mexp2int1,x)
mexp2int2 = mexp2int2 + mexp2int2.subs(x,l)
#print('EIv(x) = ', mexp2int2,'\n')

#print('Solve for C at L-a where deflection is 0') #above first azxle
#print('EIv(l-a) = 0 = ',simplify(mexp2int2.subs(x,(l-a))), '\n')

#----------------------------------Solving for Supports--------------------------
#Solution C
mexp2int2=expand(mexp2int2)
eqC = mexp2int2.subs([(C,0),(x,(l-a))])/(-1*mexp2int2.subs([(x,(l-a)),(w,0),(C,1)]))
#print('C = ',eqC,'\n')

#print('Solve for B using Sum of Forces')
#print('0 = w*l + ', eqA, ' + B ', eqC,'\n')
sumf = w*l + eqA + eqC + B
#print('0 = ', sumf,'\n')
sumf = expand(sumf)

#Solution B
eqB = simplify(sumf.subs(B,0)/sumf.subs([(w,0),(B,1)]))
#print('B = ',eqB,'\n')

#Solution A
#print('Solve for A using B')
solA = eqA.subs(B,eqB)
#print(solA,'\n')

reactionA = N(solA.subs([(w,weight_per_ft),(l,length),(a,axle_offset)]),8)
reactionC = N(eqC.subs([(w,weight_per_ft),(l,length),(a,axle_offset)]),8)
reactionB = (weight_per_ft*length) - reactionA - reactionC
reactionMa = ma.subs([(A,reactionA),(B,reactionB),(w,weight_per_ft),(l,length),(a,axle_offset)])

print('REACTIONARY FORCES')
print('A = ',reactionA,"lbs")
print('B = ',reactionB, "lbs")
print('C = ',reactionC,"lbs")
print('Ma = ',reactionMa,"lbs*ft")