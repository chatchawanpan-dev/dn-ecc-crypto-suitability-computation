# Computation Archive for the JDMSC D(n)-Triple Elliptic-Curve Screen

This repository contains the computation files and numerical outputs for the
manuscript:

**Legendre--Twist Point Counts for Split-Cubic Curves from D(n)-Triples: A
Bounded Arithmetic Screen**

The repository is intended as the reproducibility archive for the manuscript.
It contains the Magma source files, captured transcripts, CSV exports, checksum
manifest, and environment notes needed to reproduce the finite-field arithmetic
screen.

## Repository Layout

```text
.
├── Magma/
│   ├── enum_triples.m
│   ├── curve_Eabc.m
│   ├── crypto_tests.m
│   ├── all_hits.m
│   ├── worked_example.m
│   └── export_all_rows_csv.m
├── transcripts/
│   ├── enum_triples.out
│   ├── curve_Eabc.out
│   ├── crypto_tests.out
│   ├── all_hits.out
│   └── worked_example.out
├── data/
│   ├── all_rows_part01.csv
│   ├── all_rows_part02.csv
│   ├── all_rows_part03.csv
│   ├── all_rows_part04.csv
│   ├── h4_hits.csv
│   └── csv_export_summary.txt
├── ENVIRONMENT.txt
├── SHA256SUMS.txt
└── README.md
```

## What the Computation Does

The computation studies the split-cubic elliptic curves

```text
E_{a,b,c}^{(n)} : y^2 = (a x + n)(b x + n)(c x + n)
```

arising from positive integral `D(n)`-triples. The finite screen uses

```text
n in {-3, 3, 5, 8, 12, 20}
1 <= a < b < c <= 100
p_k = NextPrime(10^k), 3 <= k <= 12
```

For each sampled good reduction, the scripts compute or verify:

- the integral Weierstrass model associated with the split cubic;
- finite-field point counts;
- largest prime divisor `q` of the group order;
- cofactor `h`, relative to `q`;
- exact embedding degree `k = ord_q(p)`;
- anomalous/non-anomalous status;
- the Legendre parameter and twist factor from the manuscript;
- the Legendre trace formula used to interpret the point counts;
- all cases in the bounded scan with cofactor `h = 4`;
- the worked example for `{1, 8, 17}` with `n = 8`.

The computation is a deterministic finite arithmetic screen. It is not a
cryptographic parameter-generation procedure and is not a statistical sample
for asymptotic density claims.

## File Details

### Magma source files

| File | Purpose |
|---|---|
| `Magma/enum_triples.m` | Enumerates all positive integral `D(n)`-triples in the stated bound. |
| `Magma/curve_Eabc.m` | Builds the integral Weierstrass models for the representative triples and prints basic invariants. |
| `Magma/crypto_tests.m` | Computes representative finite-field reductions over the ten test primes. |
| `Magma/all_hits.m` | Runs the full `164 x 10` bounded scan and records all rows with cofactor `h = 4`. |
| `Magma/worked_example.m` | Verifies the worked example `{1,8,17}` with `n = 8`, including twist data and direct character-sum checks. |
| `Magma/export_all_rows_csv.m` | Exports the full row-level scan and the `h = 4` subtable as CSV files. |

### Captured transcripts

The files in `transcripts/` are the captured Magma outputs for the corresponding
scripts. The transcript-producing scripts use `SetLogFile(... : Overwrite :=
true)` so rerunning the scripts overwrites the associated transcript.

### CSV data

The four `data/all_rows_part*.csv` files contain all `1640` scanned reductions,
with the header repeated in each part. They can be reconstructed into a single
table with:

```text
awk 'FNR == 1 && NR != 1 {next} {print}' data/all_rows_part*.csv > data/all_rows.csv
```

The file `data/h4_hits.csv` contains the `45` rows for which the cofactor is
`h = 4`. The CSV columns are:

```text
n, a, b, c, p, lambda, delta, chi_delta, ap_legendre, N,
factorization_N, q, h, k, b_MOV
```

Here `lambda` and `delta` are least nonnegative representatives in `F_p`,
`N = #E(F_p)`, `factorization_N` records the prime-power factorization of `N`,
and `b_MOV = floor(k log_2 p) + 1`.

## Reproducing the Outputs

Use Magma V2.29-6 or a compatible version. From the repository root:

```text
cd Magma
magma enum_triples.m
magma curve_Eabc.m
magma crypto_tests.m
magma all_hits.m
magma worked_example.m
magma -b export_all_rows_csv.m
```

The first five commands regenerate the transcript files. The final batch command
regenerates the single CSV outputs in `data/`; the deposited row-level full
table is split into four parts for repository upload.

## Extending the Screen

To increase the enumeration bound, edit the driver variables in the relevant
Magma scripts, for example:

```text
H := 150;
Ns := [-3, 3, 5, 8, 12, 20];
```

To test larger primes, replace or extend the prime list:

```text
primes := [1009, 10007, 100003, 1000003, 10000019,
           100000007, 1000000007, 10000000019,
           100000000003, 1000000000039,
           NextPrime(2^64), NextPrime(2^96), NextPrime(2^128)];
```

Larger bounds and larger primes may require substantially longer run times.

## Checksums

Verify the deposited files with:

```text
shasum -a 256 -c SHA256SUMS.txt
```

The environment used for the submitted computations is recorded in
`ENVIRONMENT.txt`.
