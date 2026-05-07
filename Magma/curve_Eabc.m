// =====================================================================
// curve_Eabc.m
// ---------------------------------------------------------------------
// For a D(n)-triple {a,b,c}, construct the elliptic curve
//
//     E_{a,b,c}^{(n)}  :  y^2 = (a x + n)(b x + n)(c x + n)
//
// The substitution Y = (abc) y,  X = (abc) x   (which is a Q-isomorphism
// away from primes dividing abc) produces the integral model
//
//     Y^2 = X^3 + n(ab+ac+bc) X^2 + (abc) n^2 (a+b+c) X + (abc)^2 n^3.
//
// The three non-trivial rational 2-torsion points are (-abc n/a, 0) =
// (-bc n, 0), (-ac n, 0), (-ab n, 0), so the full 2-torsion subgroup
// is automatically rational:  E(Q)[2] ~ Z/2 x Z/2.
//
// Companion script for the paper
//   "Legendre--Twist Point Counts for Split-Cubic Curves from
//    D(n)-Triples: A Bounded Arithmetic Screen"
//
// Reproduction:  magma curve_Eabc.m
// =====================================================================

SetLogFile("curve_Eabc.out" : Overwrite := true);

CurveFromTriple := function(a, b, c, n)
    A := n*(a*b + a*c + b*c);
    B := a*b*c*n^2*(a + b + c);
    C := (a*b*c)^2*n^3;
    return EllipticCurve([0, A, 0, B, C]);
end function;

ReportCurve := procedure(a, b, c, n)
    printf "\n---- n = %o,  (a,b,c) = (%o,%o,%o) ----\n", n, a, b, c;
    E := CurveFromTriple(a, b, c, n);
    printf "Weierstrass coefficients [a1,a2,a3,a4,a6] = %o\n", aInvariants(E);
    printf "Discriminant                              = %o\n", Discriminant(E);
    printf "j-invariant                               = %o\n", jInvariant(E);
    T2, _ := TwoTorsionSubgroup(E);
    printf "E(Q)[2]                                   = %o\n", T2;
end procedure;

printf "=== Associated elliptic curves for representative triples ===\n";

ReportCurve( 1,  4,  7, -3);
ReportCurve( 1,  6, 13,  3);
ReportCurve( 1,  4, 11,  5);
ReportCurve( 1,  8, 17,  8);
ReportCurve( 1,  4, 13, 12);
ReportCurve( 1,  5, 16, 20);

UnsetLogFile();
quit;
