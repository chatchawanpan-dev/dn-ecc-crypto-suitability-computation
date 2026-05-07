// =====================================================================
// all_hits.m
// ---------------------------------------------------------------------
// Exhaustive h = 4 scan for the manuscript
//   "Legendre--Twist Point Counts for Split-Cubic Curves from
//    D(n)-Triples: A Bounded Arithmetic Screen"
//
// For each positive integral D(n)-triple with 1 <= a < b < c <= H and
// n in {-3,3,5,8,12,20}, reduce E_{a,b,c}^{(n)} at ten decimal-size
// test primes
// and print every row for which the resulting cofactor h equals 4.
// This script reproduces the exhaustive hit table in the paper.
//
// Reproduction: magma all_hits.m
// =====================================================================

SetLogFile("all_hits.out" : Overwrite := true);

IsIntegerSquare := function(m)
    return m ge 0 and IsSquare(m);
end function;

EnumerateDnTriples := function(n, H)
    T := [];
    for a in [1..H] do
        for b in [(a+1)..H] do
            if not IsIntegerSquare(a*b + n) then continue; end if;
            for c in [(b+1)..H] do
                if not IsIntegerSquare(a*c + n) then continue; end if;
                if not IsIntegerSquare(b*c + n) then continue; end if;
                Append(~T, [a,b,c]);
            end for;
        end for;
    end for;
    return T;
end function;

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

H := 100;
Ns := [-3, 3, 5, 8, 12, 20];
primes := [1009, 10007, 100003, 1000003, 10000019,
           100000007, 1000000007, 10000000019,
           100000000003, 1000000000039];

row_count := 0;
hit_count := 0;

printf "=== Exhaustive h = 4 scan up to H = %o ===\n", H;
printf "Columns: n, triple, p, #E(F_p), h, q, bits, k\n\n";

for n in Ns do
    T := EnumerateDnTriples(n, H);
    for abc in T do
        a := abc[1];  b := abc[2];  c := abc[3];
        E := CurveFromTriple(a, b, c, n);
        for p in primes do
            row_count +:= 1;
            AssertGoodReduction(E, p);
            Ep := ChangeRing(E, GF(p));
            N := #Ep;
            assert N mod 4 eq 0;
            assert N ne p;
            lambda, delta, chiDelta, apL := LegendreTraceData(a, b, c, n, p, N);
            q, fac := LargestPrimeFactor(N);
            assert IsPrime(q);
            h := N div q;
            assert h eq N div q;
            k := EmbeddingDegreeExact(q, p);
            assert k eq Order(Integers(q)!(p mod q));
            if h eq 4 then
                hit_count +:= 1;
                assert N eq 4*q;
                assert IsPrime(q);
                printf "%3o  (%o,%o,%o)  %8o  %8o  %3o  %8o  %2o  %o\n",
                       n, a, b, c, p, N, h, q, Ilog2(q)+1, k;
            end if;
        end for;
    end for;
end for;

printf "\nTotal rows scanned = %o\n", row_count;
printf "Total h = 4 hits   = %o\n", hit_count;

UnsetLogFile();
quit;
