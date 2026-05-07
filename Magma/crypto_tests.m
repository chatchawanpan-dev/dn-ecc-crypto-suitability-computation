// =====================================================================
// crypto_tests.m
// ---------------------------------------------------------------------
// For each representative triple, reduce E_{a,b,c}^{(n)} modulo a ladder
// of test primes p and report:
//   * #E(F_p),                              (point count)
//   * factorisation of #E(F_p),             (cofactor h, largest prime q)
//   * whether |E(F_p)| = p                   (anomalous ? Smart attack)
//   * embedding degree k = ord_q(p)          (MOV / FR)
//
// The full 2-torsion over Q forces 4 | #E(F_p) for every odd good prime p,
// so these reductions are automatically non-anomalous. For the test primes
// used here, the resulting cofactor is necessarily at least 4.
//
// Companion script for the paper
//   "Legendre--Twist Point Counts for Split-Cubic Curves from
//    D(n)-Triples: A Bounded Arithmetic Screen"
//
// Reproduction:  magma crypto_tests.m
// =====================================================================

SetLogFile("crypto_tests.out" : Overwrite := true);

CurveFromTriple := function(a, b, c, n)
    A := n*(a*b + a*c + b*c);
    B := a*b*c*n^2*(a + b + c);
    C := (a*b*c)^2*n^3;
    return EllipticCurve([0, A, 0, B, C]);
end function;

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

CryptoReport := procedure(a, b, c, n, E, p)
    AssertGoodReduction(E, p);
    Ep := ChangeRing(E, GF(p));
    N  := #Ep;
    assert N mod 4 eq 0;
    assert N ne p;
    lambda, delta, chiDelta, apL := LegendreTraceData(a, b, c, n, p, N);

    printf "  p = %8o  #E(F_p) = %-12o  =  %o\n",
           p, N, Factorisation(N);
    q, fac := LargestPrimeFactor(N);
    assert IsPrime(q);
    h := N div q;
    assert h eq N div q;
    k := EmbeddingDegreeExact(q, p);
    assert k eq Order(Integers(q)!(p mod q));
    anomalous := (N eq p);
    printf "                        q = %o (%o-bit),  h = %o,  k = %o,  anomalous = %o\n",
           q, Ilog2(q)+1, h, k, anomalous;
    printf "                        Legendre check: chi(delta) = %o,  a_p(L_lambda) = %o\n",
           chiDelta, apL;
end procedure;

CryptoBatch := procedure(a, b, c, n, primes)
    printf "\n==== n = %o,  (a,b,c) = (%o,%o,%o) ====\n", n, a, b, c;
    E := CurveFromTriple(a, b, c, n);
    for p in primes do
        CryptoReport(a, b, c, n, E, p);
    end for;
end procedure;

primes := [1009, 10007, 100003, 1000003, 10000019,
           100000007, 1000000007, 10000000019,
           100000000003, 1000000000039];

CryptoBatch( 1,  4,  7, -3, primes);
CryptoBatch( 1,  6, 13,  3, primes);
CryptoBatch( 1,  4, 11,  5, primes);
CryptoBatch( 1,  8, 17,  8, primes);
CryptoBatch( 1,  4, 13, 12, primes);
CryptoBatch( 1,  5, 16, 20, primes);

UnsetLogFile();
quit;
