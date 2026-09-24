# VALD parser memory audit

## Scope

This audit targets the temporary memory peak while loading the canonical
v1.2.0 cool, metal-rich 800 A line list.  The input is a counted VALD
`long + extract_stellar` file containing 1,258,596 spectral lines and occupying
0.535 GiB on disk.

The implementation adds a deliberately narrow streaming fast path for this
unambiguous format.  Instead of materializing the whole file, its line-record
slices, and a second joined string at once, the reader parses 20,000 spectral
lines per chunk and concatenates the resulting data frames.  Reference keys
are accumulated per chunk.  Short-format, `extract_all`, and other VALD
variants continue to use the unchanged legacy parser.

## Result

Fresh-process 10 ms RSS sampling measured:

| Metric | legacy parser | streaming parser | Difference |
|---|---:|---:|---:|
| Peak RSS | 4,186,488,832 B | 1,856,094,208 B | -2,330,394,624 B (-55.7%) |
| Peak RSS (GiB) | 3.899 | 1.729 | -2.170 |
| Parse wall time | not isolated in the memory run | 6.36 s | -- |

The final data frame is unchanged in size (231,581,796 B shallow;
873,257,734 B deep), contains all 1,258,596 lines, and spans the expected
4650.0003--5749.9987 A range.

## Verification and stop condition

A forced multi-chunk regression test compares the complete frame, dtypes,
atmosphere, abundances, references, units, and format metadata exactly against
the legacy parser on the repository long-form fixture.  All 11 VALD tests
pass.

This removes more than half of the parser peak without changing the parser for
ambiguous or less common VALD formats.  Further reductions would require
changing object-column storage or preallocating the final frame, both of which
have a larger compatibility cost and are outside this release's stop
condition.  A full 800 A synthesis was not rerun for this parser-only change
to avoid an unnecessary battery-intensive benchmark; end-to-end synthesis
remains covered by the existing line-state and release-candidate measurements.
