// =====================================================================
// export_all_rows_csv.m
// ---------------------------------------------------------------------
// Generate machine-readable CSV files for the bounded exhaustive scan:
//   ../all_rows.csv  -- all 1640 scanned reductions
//   ../h4_hits.csv   -- the subtable with h = 4
//
// Reproduction from this directory:
//   magma -b export_all_rows_csv.m
// =====================================================================

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

FactorisationString := function(fac)
    if #fac eq 0 then
        return "1";
    end if;

    s := "";
    for i in [1..#fac] do
        term := Sprintf("%o^%o", fac[i][1], fac[i][2]);
        if i eq 1 then
            s := term;
        else
            s := s cat ";" cat term;
        end if;
    end for;
    return s;
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

    return Integers()!lambda, Integers()!delta, chiDelta, apL;
end function;

TransferFieldBits := function(p, k)
    R := RealField(200);
    return Integers()!(Floor((R!k) * Log(R!p) / Log(R!2)) + 1);
end function;

H := 100;
Ns := [-3, 3, 5, 8, 12, 20];
primes := [1009, 10007, 100003, 1000003, 10000019,
           100000007, 1000000007, 10000000019,
           100000000003, 1000000000039];

AllRows := Open("../all_rows.csv", "w");
H4Rows := Open("../h4_hits.csv", "w");
Summary := Open("../csv_export_summary.txt", "w");

header := "n,a,b,c,p,lambda,delta,chi_delta,ap_L_lambda,N,factorization_N,q,h,k,b_mov\n";
fprintf AllRows, header;
fprintf H4Rows, header;

row_count := 0;
hit_count := 0;

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
            facString := FactorisationString(fac);
            assert IsPrime(q);
            h := N div q;
            k := EmbeddingDegreeExact(q, p);
            bMOV := TransferFieldBits(p, k);

            fprintf AllRows, "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
                    n, a, b, c, p, lambda, delta, chiDelta, apL, N, facString, q, h, k, bMOV;
            if h eq 4 then
                hit_count +:= 1;
                fprintf H4Rows, "%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o,%o\n",
                        n, a, b, c, p, lambda, delta, chiDelta, apL, N, facString, q, h, k, bMOV;
            end if;
        end for;
    end for;
end for;

fprintf Summary, "Total rows scanned = %o\n", row_count;
fprintf Summary, "Total h = 4 hits   = %o\n", hit_count;

quit;
