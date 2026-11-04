restart:
with(LinearAlgebra):
Digits := 15:

# Independent direct integral check of check_effective_potential.
# Set parameters to the values printed by the Fortran test.
nh := 5:
eta := 1.3544: rho := 1.: sigma := -1.: gam := 0.5:
kappa1 := 1.: kappa2 := 1.: bath12 := 0.1:
q := Vector([1.,2.,3.]):
B := Matrix([[1.,0.1,0.15],[0.1,2.,0.15],
             [0.15,0.15,3.]]):
W := Matrix([[kappa1/2,bath12/2],
             [bath12/2,kappa2/2]]):
c := Vector(nh,1):

# The Fortran Hermite basis uses physicists' HermiteH,
# argument sqrt(B[1,1])*(x-q[1]), normalized by sqrt(2^n*n!).
H := proc(i,z)
    global B, q;
    return HermiteH(i-1,sqrt(B[1,1])*(z-q[1]))
           /sqrt(2^(i-1)*factorial(i-1));
end proc:

# Integrate the two bath coordinates exactly, then integrate x.
A := SubMatrix(B,2..3,2..3):
av := Vector([B[2,1],B[3,1]]):
Ai := MatrixInverse(A):
alpha := simplify(B[1,1]-DotProduct(av,Ai.av)):
Avec := Ai.av:
mu := Vector([seq(q[k+1]-Avec[k]*(x-q[1]),k=1..2)]):
N2 := sqrt(Determinant(B)/Pi^3):
Y0 := Pi/sqrt(Determinant(A)):
weight := N2*Y0*exp(-alpha*(x-q[1])^2):

# The potential is V = Vx + gam*x*y1 + y^T W y.
# Conditional Gaussian y|x has covariance Ai/2.
Vybar := expand(DotProduct(mu,W.mu)+Trace(W.Ai)/2):
Vbar := expand(x^4/(16*eta)+rho*x^3+sigma*x^2/2
               +gam*x*mu[1]+Vybar):
Dbar := Vector([
  x^3/(4*eta)+3*rho*x^2+sigma*x+gam*mu[1],
  gam*x+2*(W.mu)[1],
  2*(W.mu)[2]
]):
DDbar := Matrix([
  [3*x^2/(4*eta)+6*rho*x+sigma,gam,0],
  [gam,2*W[1,1],2*W[1,2]],
  [0,2*W[2,1],2*W[2,2]]
]):

# Polynomial times Gaussian: Maple evaluates these as direct integrals.
IntegrateX := proc(f)
    global weight, x;
    return evalf(int(expand(f)*weight,x=-infinity..infinity));
end proc:

S := Matrix(nh,nh): V0 := Matrix(nh,nh):
for i from 1 to nh do
  for j from 1 to nh do
    S[i,j] := IntegrateX(H(i,x)*H(j,x)):
    V0[i,j] := IntegrateX(H(i,x)*H(j,x)*Vbar):
  end do:
end do:

V1 := Vector(3): V2 := Matrix(3,3):
for k from 1 to 3 do
  V1[k] := add(add(c[i]*c[j]*IntegrateX(
               H(i,x)*H(j,x)*Dbar[k]),j=1..nh),i=1..nh):
  for l from 1 to 3 do
    V2[k,l] := add(add(c[i]*c[j]*IntegrateX(
                 H(i,x)*H(j,x)*DDbar[k,l]),j=1..nh),i=1..nh):
  end do:
end do:

printf("S00 (overlap matrix):\n"):
print(evalf(S)):
printf("Norm c^T S c = %g\n",evalf(DotProduct(c,S.c))):
printf("V0 matrix:\n"):
print(evalf(V0)):
printf("V1 vector:\n"):
print(evalf(V1)):
printf("V2 matrix:\n"):
print(evalf(V2)):
printf("c^T V0 c = %g\n",evalf(DotProduct(c,V0.c))):

# Isolate the cubic contribution without relying on Fortran formulas.
X1 := Matrix(nh,nh): X2 := Matrix(nh,nh): X3 := Matrix(nh,nh):
for i from 1 to nh do
  for j from 1 to nh do
    X1[i,j] := IntegrateX(H(i,x)*H(j,x)*x):
    X2[i,j] := IntegrateX(H(i,x)*H(j,x)*x^2):
    X3[i,j] := IntegrateX(H(i,x)*H(j,x)*x^3):
  end do:
end do:
printf("cubic delta V0:\n"): print(evalf(rho*X3)):
printf("cubic delta V1_X = %g\n",
       evalf(3*rho*DotProduct(c,X2.c))):
printf("cubic delta V2_XX = %g\n",
       evalf(6*rho*DotProduct(c,X1.c))):
