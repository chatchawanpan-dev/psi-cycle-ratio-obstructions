// magma_verify_and_explore.m
// Independent Magma verification and exploration for
//   psi_{d,a,b}(z) = z/(a z^d + b).
//
// Run from the repository root with:
//   magma -b code/magma_verify_and_explore.m

Q := Rationals();

print "============================================================";
print "Magma verification for psi_{d,a,b}(z)=z/(a z^d+b)";
print "============================================================";

// -------------------------------------------------------------------------
// 1. Period-2 identity for d=3.
// -------------------------------------------------------------------------
Puv<u,v> := PolynomialRing(Q, 2);
Phi2_uv := v*u^2 + (2*v^2 + v - 1)*u + v^2*(v + 1);
assert (u+v)^2 - u - v*(u+v)^3 eq -(u+v-1)*Phi2_uv;
print "PASS: period-2 quotient identity in Z[u,v] for d=3.";

K2<A,B> := FunctionField(Q, 2);
Rw<W> := PolynomialRing(K2);
Phi2_w := A^2*B*W^2 + A*(2*B^2 + B - 1)*W + B^2*(B + 1);
assert Discriminant(Phi2_w) eq -A^2*(3*B - 1)*(B + 1);
print "PASS: discriminant(Phi_2^*(W)) = -A^2(3B-1)(B+1).";

// -------------------------------------------------------------------------
// 2. Rational period-2 parameter locus for d=3.
// -------------------------------------------------------------------------
Krt<R,Tpar> := FunctionField(Q, 2);
Dpar := R^2 + R + 1;
Apar := (R + 1)/(R*Dpar*Tpar^3);
Bpar := R/Dpar;

function PsiParam(X)
    return X/(Apar*X^3 + Bpar);
end function;

assert PsiParam(Tpar) eq R*Tpar;
assert PsiParam(R*Tpar) eq Tpar;
print "PASS: rational period-2 parameterization t <-> r*t verified.";

KR<Rc> := FunctionField(Q);
M := Matrix(KR, 2, 2, [1, 1, Rc^3, 1]);
assert Determinant(M) eq 1 - Rc^3;
Aconv := (1/Rc - Rc)/(1 - Rc^3);
Bconv := (Rc - Rc^2)/(1 - Rc^3);
assert Aconv eq (Rc + 1)/(Rc*(Rc^2 + Rc + 1));
assert Bconv eq Rc/(Rc^2 + Rc + 1);
assert Aconv + Bconv eq 1/Rc;
assert Rc^3*Aconv + Bconv eq Rc;
print "PASS: converse linear solve gives A=(r+1)/(r(r^2+r+1)), b=r/(r^2+r+1).";

// -------------------------------------------------------------------------
// 3. General period-2 cycle-ratio identities.
// -------------------------------------------------------------------------
Kt2<Tau> := FunctionField(Q);
for d in [2..8] do
    Btau := (Tau^(d-1) - Tau)/(Tau^d - 1);
    assert Btau*Tau^d - Tau^(d-1) + Tau - Btau eq 0;
    assert (Tau^d - 1)^2 - (Tau^(d-1) - Tau)^2
        eq (Tau^2 - 1)*(Tau^(2*d-2) - 1);
end for;
print "PASS: general period-2 cycle-ratio identities checked for d=2..8.";

Prs<rr,ss> := PolynomialRing(Q, 2);
E3 := rr^4*ss^4 - rr^3*ss^2 - rr^2*ss^3 + rr^2 - rr*ss + ss^2;
Elim3 := (rr^3*ss-rr)*(rr^3*ss^3-1)-(rr^2*ss^2-rr)*(rr^3-1);
assert Elim3 eq rr^2*E3;
assert 4*E3 eq (2*rr^2*ss^2-rr-ss)^2 + 3*(rr-ss)^2;
print "PASS: period-3 ratio elimination and sum-of-squares identity verified.";

// -------------------------------------------------------------------------
// 4. Leading/constant coefficient check for the odd-d period-2 reduction.
//    The manuscript proves the general formula; this checks several odd d.
// -------------------------------------------------------------------------
K1<B1> := FunctionField(Q);
Ru<U> := PolynomialRing(K1);
for d in [3,5,7,9,11] do
    E := (U+B1)^(d-1) - U - B1*(U+B1)^d;
    q, r := Quotrem(E, -(U+B1-1));
    assert r eq 0;
    assert LeadingCoefficient(q) eq B1;
    assert Coefficient(q, 0) eq B1^(d-1)*(B1+1);
end for;
print "PASS: odd-d coefficient checks for d = 3,5,7,9,11.";

// -------------------------------------------------------------------------
// 5. Period-3 dynatomic factorization for d=3.
// -------------------------------------------------------------------------
K3<A3,B3> := FunctionField(Q, 2);
Rz<Z> := PolynomialRing(K3);
Fz := FieldOfFractions(Rz);
psi := Fz!Z / (A3*(Fz!Z)^3 + B3);
psi2 := Evaluate(psi, psi);
psi3 := Evaluate(psi, psi2);
num3 := Numerator(psi3 - Fz!Z);
P1 := -Z*(A3*Z^3 + B3 - 1);
quot3, rem3 := Quotrem(num3, P1);
assert rem3 eq 0;

Q1_Z := A3^2*Z^6 + A3*(2*B3+1)*Z^3 + (B3^2+B3+1);
WZ := Z^3;
S_Z :=
    A3^6*B3^4*WZ^6
    + A3^5*(6*B3^5-B3^2)*WZ^5
    + A3^4*(15*B3^6-2*B3^3+1)*WZ^4
    + A3^3*(20*B3^7+B3)*WZ^3
    + A3^2*(15*B3^8+2*B3^5+B3^2)*WZ^2
    + A3*(6*B3^9+B3^6)*WZ
    + B3^10;

// Magma normalizes the numerator over Q(A3,B3), so compare up to the
// normalizing scalar A3^9 B3^4.
assert A3^9*B3^4*quot3 eq Q1_Z*S_Z;
print "PASS: period-3 factorization Q1(Z^3)*S(Z^3) verified.";

Kdisc<A4,B4> := FunctionField(Q, 2);
RQ<W4> := PolynomialRing(Kdisc);
Q1_W := A4^2*W4^2 + A4*(2*B4+1)*W4 + (B4^2+B4+1);
assert Discriminant(Q1_W) eq -3*A4^2;
print "PASS: discriminant(Q1) = -3 A^2.";

// -------------------------------------------------------------------------
// 6. Rational-b exploration for period 2.
// -------------------------------------------------------------------------
Kt<T> := FunctionField(Q);
b_param := -2*(T+1)/(T^2+3);
y_param := 1 + T*b_param;
assert y_param^2 eq -(3*b_param-1)*(b_param+1);
print "PASS: discriminant-square conic parametrized by";
print "      b = -2(t+1)/(t^2+3), y = 1 + t b.";

function Psi3(x, a, b)
    return x/(a*x^3 + b);
end function;

a_ex := Q!1/6;
b_ex := Q!(-2/3);
assert Psi3(Q!1, a_ex, b_ex) eq Q!(-2);
assert Psi3(Q!(-2), a_ex, b_ex) eq Q!1;
print "EXAMPLE: for a=1/6, b=-2/3, the rational 2-cycle is 1 <-> -2.";

function IsRatCube(r)
    n := Numerator(r);
    d := Denominator(r);
    okn, rn := IsPower(n, 3);
    if not okn then
        return false, Q!0;
    end if;
    okd, rd := IsPower(d, 3);
    if not okd then
        return false, Q!0;
    end if;
    return true, Q!rn/Q!rd;
end function;

Pq<Wq> := PolynomialRing(Q);
found_a1 := [];
height_bound := 40;
for den in [1..height_bound] do
    for num in [-height_bound..height_bound] do
        if num ne 0 and GCD(Abs(num), den) eq 1 then
            bv := Q!num/den;
            phi := bv*Wq^2 + (2*bv^2 + bv - 1)*Wq + bv^2*(bv+1);
            for rt in Roots(phi) do
                ok, z0 := IsRatCube(rt[1]);
                if ok and z0 ne 0 then
                    Append(~found_a1, <bv, rt[1], z0>);
                end if;
            end for;
        end if;
    end for;
end for;
print "SEARCH: a=1, rational b of height <= 40 with rational-cube period-2 root:";
print found_a1;

// -------------------------------------------------------------------------
// 7. Search the displayed period-3 sextic for rational roots at a=1.
// -------------------------------------------------------------------------
function SexticS(bv)
    return
        bv^4*Wq^6
        + (6*bv^5-bv^2)*Wq^5
        + (15*bv^6-2*bv^3+1)*Wq^4
        + (20*bv^7+bv)*Wq^3
        + (15*bv^8+2*bv^5+bv^2)*Wq^2
        + (6*bv^9+bv^6)*Wq
        + bv^10;
end function;

bad_sextic_roots := [];
bad_sextic_cubes := [];
for bv_int in [-100..100] do
    if bv_int ne 0 then
        bv := Q!bv_int;
        roots := Roots(SexticS(bv));
        if #roots gt 0 then
            Append(~bad_sextic_roots, <bv_int, roots>);
        end if;
        for rt in roots do
            ok, z0 := IsRatCube(rt[1]);
            if ok and z0 ne 0 then
                Append(~bad_sextic_cubes, <bv_int, rt[1], z0>);
            end if;
        end for;
    end if;
end for;
print "SEARCH: rational roots of S_{1,b}(W) for integer b in [-100,100]\\{0}:";
print bad_sextic_roots;
print "SEARCH: rational-cube roots of S_{1,b}(W) in the same range:";
print bad_sextic_cubes;

print "============================================================";
print "Magma verification complete.";
print "============================================================";

quit;
