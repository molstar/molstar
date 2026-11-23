# Drift: locations-and-queries.md

## 1. The row form and a pair form

**What**: the design's first forms list `(row :id k)` giving Keys ("with alignments") and, among Keys, `(pair …)`. The
code has `(row :key k)` (and a list of keys) under an alignment, giving a Location on its rows system (a side of a
rect), with `seq.query.row.unknown` and `.pending`; there is no `pair` form: a pairs annotation answers `items` (its
pairs as Keys), as the design's own later decision ("Pairs answer `items`") says.

**Where**: `packages/seq/core/src/query/` (the `row` resolver); the core spec, Queries, Forms and `row`.

**Options**: (a) the design follows the code: rows are a system (EC-35), so a row form is a place on it, and pairs are
items; (b) add a Keys row form and a pair form.

**Recommendation**: (a). A row as a position is what rects, selection and the rows layer read; nothing asks for row
keys.

## 2. How an extension's form prints

**What**: the design says an extension form's params print "in declared order". The code prints every form by the
address grammar, keys sorted (core's printer), whatever the declaration says (core spec, Queries, Extensions).

**Where**: `@molstar/core/repr/tree` (`printAddress`); `packages/seq/core/src/query/extensions.ts`.

**Options**: (a) the design follows the code: one canonical spelling for every form, so equal addresses print alike;
(b) a printer that reads the extension's declaration.

**Recommendation**: (a): an address's identity as text must not depend on whether the reader has the extension.

**Resolved** (the author, 2026-11-23, as recommended): both items (a): the design's form table and the printing rule now say the code's way (`(row :key k)` a Location on the rows system, pairs answer `items`; every form printed sorted by key); the two rows point at core › Forms and core › Extensions.
