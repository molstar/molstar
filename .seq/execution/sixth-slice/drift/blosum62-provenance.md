# Drift: what the `blosum62` scheme tells a picker

**What**: `theme/provenance.ts` describes `blosum62` to a picker as "Jalview's BLOSUM62 score against the column's
consensus". The scheme itself (`theme/schemes.ts`) is a rules form over the `matches-reference` measure: a letter
equal to the reference row's at that column (1 or 0, NaN without a reference). No BLOSUM62 matrix is read and no
consensus is computed for it. The core spec (Themes, "Column schemes") says what the scheme does. themes.md lists the
scheme by name only.

**Where**: `packages/seq/core/src/theme/provenance.ts` (the `blosum62` entry); `packages/seq/core/src/theme/schemes.ts`
(the scheme's form); core spec, Themes, Column schemes.

**Options**: (a) correct the provenance text to say what is drawn ("matches the reference row's letter; no
substitution score"), and perhaps rename the scheme's label; (b) build the scheme as named, with a BLOSUM62 score of
each letter against the column's consensus (a matrix in core's data, a `consensus` measure on the rows layer).

**Recommendation**: (a) now, as a fix lane's item: a picker must never claim a score the picture does not show. (b)
when an alignment view asks for substitution scores; then the name is earned and the spec's row changes with it.
