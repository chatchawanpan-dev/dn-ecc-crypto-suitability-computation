// =====================================================================
// enum_triples.m
// ---------------------------------------------------------------------
// Enumerate D(n)-triples {a,b,c} with 1 <= a < b < c <= H, i.e. ordered
// triples of positive integers with the property that ab+n, ac+n, bc+n
// are all perfect squares.
//
// Companion script for the paper
//   "Legendre--Twist Point Counts for Split-Cubic Curves from
//    D(n)-Triples: A Bounded Arithmetic Screen"
//
// Reproduction:  magma enum_triples.m
// Expected runtime on H=100 : well under 1 second.
// =====================================================================

SetLogFile("enum_triples.out" : Overwrite := true);

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

// ---------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------
printf "=== D(n)-triple enumeration up to H = 100 ===\n";

H   := 100;
Ns  := [-3, 3, 5, 8, 12, 20];

for n in Ns do
    T := EnumerateDnTriples(n, H);
    printf "n = %4o  :  %o D(n)-triples\n", n, #T;
    if #T gt 0 then
        printf "  first triple  : %o\n", T[1];
        printf "  last  triple  : %o\n", T[#T];
    end if;
end for;

// Show the full list for n = -3 as an illustration.
printf "\nFull list for n = -3:\n";
T := EnumerateDnTriples(-3, H);
for i in [1..#T] do
    printf "  %2o. %o\n", i, T[i];
end for;

UnsetLogFile();
quit;
