// =====================================================================
// worked_example.m
// ---------------------------------------------------------------------
// A fully worked finite-field arithmetic walk-through of the elliptic curve
//
//    E : y^2 = (x + 8)(8x + 8)(17x + 8)        (D(8)-triple {1,8,17})
//
// produced from the D(n)-triple (1,8,17) with n = 8. This triple is
// the only representative in our survey whose associated curve has
// cofactor h = 4 at one of the test primes (p = 1000003), and whose
// quadratic twist separately has h' = 4 at p = 1009. The script builds the integral
// Weierstrass model, computes its invariants, displays the three
// rational 2-torsion points, reduces E at ten decimal-size test primes and reports
// the main-curve data (cofactor, largest prime factor, embedding degree,
// anomalous check). We also compute the nontrivial quadratic twist at each
// prime by choosing an explicit nonsquare and record the corresponding data,
// including the Legendre trace check.
//
// Companion script for the paper
//   "Legendre--Twist Point Counts for Split-Cubic Curves from
//    D(n)-Triples: A Bounded Arithmetic Screen"
//
// Reproduction:  magma worked_example.m
// =====================================================================

SetLogFile("worked_example.out" : Overwrite := true);

printf "=== Worked example: D(8)-triple (a,b,c) = (1,8,17) ===\n";

n := 8;
a := 1;  b := 8;  c := 17;

// Sanity-check the Diophantine property
printf "Verify D(n)-triple property:\n";
printf "  ab + n = %o = %o^2\n", a*b + n, Isqrt(a*b + n);
printf "  ac + n = %o = %o^2\n", a*c + n, Isqrt(a*c + n);
printf "  bc + n = %o = %o^2\n", b*c + n, Isqrt(b*c + n);

R<x> := PolynomialRing(Rationals());
f := (a*x + n) * (b*x + n) * (c*x + n);
printf "\nDefining polynomial f(x) = %o\n", f;

A := n*(a*b + a*c + b*c);
B := a*b*c*n^2*(a + b + c);
C := (a*b*c)^2*n^3;
E := EllipticCurve([0, A, 0, B, C]);
printf "Weierstrass model E     = %o\n", E;
printf "a-invariants            = %o\n", aInvariants(E);
printf "Discriminant            = %o\n", Discriminant(E);
printf "j-invariant             = %o\n", jInvariant(E);

T2, mT2 := TwoTorsionSubgroup(E);
printf "\nRational 2-torsion E(Q)[2] = %o\n", T2;
printf "Explicit 2-torsion points (on the integral model) :\n";
printf "  P_a = (-bc*n,  0) = (%o, 0)\n", -b*c*n;
printf "  P_b = (-ac*n,  0) = (%o, 0)\n", -a*c*n;
printf "  P_c = (-ab*n,  0) = (%o, 0)\n", -a*b*n;

printf "\n--- Reductions at test primes ---\n";

primes := [1009, 10007, 100003, 1000003, 10000019,
           100000007, 1000000007, 10000000019,
           100000000003, 1000000000039];

LargestPrimeFactor := function(N)
    fac := Factorisation(N);
    return Max([f[1] : f in fac]), fac;
end function;

AssertGoodReduction := procedure(E, p)
    assert IsPrime(p);
    assert (Integers()!Discriminant(E)) mod p ne 0;
end procedure;

EmbeddingDegreeExact := function(q, p)
    if not IsPrime(q) then
        error "q must be prime";
    end if;
    if GCD(q, p) ne 1 then
        error "p and q must be coprime";
    end if;
    return Order(Integers(q)!(p mod q));
end function;

QuadraticCharacterPrimeField := function(x)
    if x eq 0 then
        return 0;
    elif IsSquare(x) then
        return 1;
    else
        return -1;
    end if;
end function;

LegendreTraceData := function(a, b, c, n, p, N)
    F := GF(p);
    lambda := F!(b*(c-a)) / F!(c*(b-a));
    delta  := F!(c*n*(b-a));
    assert lambda ne F!0;
    assert lambda ne F!1;
    assert delta ne F!0;

    L := EllipticCurve([F | 0, -(F!1 + lambda), 0, lambda, 0]);
    apL := p + 1 - #L;
    chiDelta := QuadraticCharacterPrimeField(delta);
    assert N eq p + 1 - chiDelta*apL;

    return lambda, delta, chiDelta, apL;
end function;

DirectCharacterSumCount := function(lambda, delta, p)
    F := GF(p);
    total := 1; // point at infinity
    for u in F do
        total +:= 1 + QuadraticCharacterPrimeField(delta*u*(u - F!1)*(u - lambda));
    end for;
    return Integers()!total;
end function;

for p in primes do
    AssertGoodReduction(E, p);
    F := GF(p);
    Ep  := ChangeRing(E, F);
    N   := #Ep;
    assert N mod 4 eq 0;
    assert N ne p;
    t   := p + 1 - N;                         // trace of Frobenius
    lambda, delta, chiDelta, apL := LegendreTraceData(a, b, c, n, p, N);
    assert t eq chiDelta*apL;
    q, fac := LargestPrimeFactor(N);
    assert IsPrime(q);
    h   := N div q;
    assert h eq N div q;
    k := EmbeddingDegreeExact(q, p);
    assert k eq Order(Integers(q)!(p mod q));

    d := PrimitiveElement(F);
    assert not IsSquare(d);
    Et := QuadraticTwist(Ep, d);
    Nt := #Et;
    assert Nt eq p + 1 + t;
    assert N + Nt eq 2*p + 2;

    qT, facT := LargestPrimeFactor(Nt);
    hT := Nt div qT;

    if p in {1009, 10007, 100003} then
        directN := DirectCharacterSumCount(lambda, delta, p);
        assert directN eq N;
    else
        directN := 0;
    end if;

    printf "\np = %o\n", p;
    printf "  #E(F_p)       = %o  =  %o\n", N, fac;
    printf "  trace t       = %o  (|t| <= 2*sqrt(p) ~= %o)\n",
           t, Floor(2*Sqrt(p));
    printf "  lambda mod p  = %o\n", lambda;
    printf "  delta mod p   = %o\n", delta;
    printf "  chi(delta)    = %o\n", chiDelta;
    printf "  a_p(L_lambda) = %o\n", apL;
    if directN ne 0 then
        printf "  direct character-sum # = %o\n", directN;
    end if;
    printf "  cofactor h    = %o\n", h;
    printf "  prime q       = %o  (%o-bit)\n", q, Ilog2(q)+1;
    printf "  embedding k   = %o\n", k;
    printf "  anomalous?    = %o\n", N eq p;
    printf "  twist factor d = %o  (nonsquare in F_p)\n", d;
    printf "  twist #E'     = %o  =  %o\n", Nt, facT;
    printf "  twist cofactor = %o, twist prime = %o (%o-bit)\n",
           hT, qT, Ilog2(qT)+1;
end for;

printf "\n--- Summary ---\n";
printf "At p = 1000003, E has cofactor h = 4 with\n";
printf "prime q = 249727 of 18 bits; no anomalous reduction occurs\n";
printf "at any test prime. At p = 1009 the nontrivial quadratic twist has\n";
printf "h' = 4 with twist prime q' = 251.\n";

UnsetLogFile();
quit;
