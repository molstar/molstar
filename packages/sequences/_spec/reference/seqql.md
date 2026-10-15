# SeqQL reference

## Selections

### all
`seq.generator.all` · generator · result: selection
`(all)`: every unit and alignment row, one group each.

Examples:
- `(and (all) (letters W))`: every tryptophan (case1)

### empty
`seq.generator.empty` · generator · result: selection
`(empty)`: nothing.

Examples:
- `(empty)`: no groups

### sel
`seq.macro.sel` · macro · result: selection
`(sel S)`: the shorthand.

| Argument | Type | Default | |
|---|---|---|---|
| `S` (positional) | text | required | refs with label ranges [10-20], index slices @9:20, item ids [#id], inputs $name; a, b a union |

Examples:
- `(sel "1brs/A[10-20]")`: residues 10-20 of 1brs/A (case7)
- `(sel "p53@0:50")`: the first 50 positions, as indices (case1)

Hints:
- `seqql.arg.missing`: sel takes the shorthand as a string: (sel "1brs/A[10-20]")
- `seqql.type`: sel takes the shorthand as a string: (sel "1brs/A[10-20]")
- `seqql.sel.syntax`: labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")
- `seqql.sel.mixed`: a bracket holds labels or #ids, not both

### systems
`seq.generator.systems` · generator · result: selection · `:where` binds: system
`(systems :of R :kind K :where T)`: whole sequences of kind unit, molecule, row or columns, one group each.

| Argument | Type | Default | |
|---|---|---|---|
| `:of` | ref | none | a molecule or complex for units, an alignment for rows; default every system |
| `:kind` | `unit` \| `molecule` \| `row` \| `columns` | `unit` | the kind of system |
| `:where` | test of a system | none | the systems kept |

Examples:
- `(systems :of 1brs)`: the units of complex 1brs, one group each (case7)

Hints:
- `seqql.arg.unknown`: it takes :of, :kind and :where; narrow positions with (and Q ...)

### rows
`seq.generator.rows` · generator · result: selection · `:where` binds: system
`(rows :of A :where T)`: an alignment's rows, one group each.

| Argument | Type | Default | |
|---|---|---|---|
| `:of` | ref | none | the alignment; default every alignment |
| `:where` | test of a system | none | the rows kept |

Examples:
- `(rows :of msa :where (>= (field identity :of traits) 0.9))`: rows at least 90% identical

Hints:
- `seqql.ref.kind`: :of takes an alignment (msa); a row alone is msa/3

### items
`seq.generator.items` · generator · result: selection · `:where` binds: item
`(items :of A :where T)`: items of annotation A (all annotations if no :of), one group each.

| Argument | Type | Default | |
|---|---|---|---|
| `:of` | ref | none | the annotation; default every annotation with items |
| `:where` | test of an item | none | the items kept |

Examples:
- `(items :of regions)`: p53's regions, one group each (case1)

Hints:
- `seqql.arg.unknown`: it takes :of and :where; items touching a selection: (pick (items :of A) :where (overlaps Q))
- `seqql.ref.kind`: values and matrices have no items: select by value with (values :of A :above x)

### values
`seq.generator.values` · generator · result: selection
`(values :of A :min x :max x :above x :below x)`: positions where values annotation A passes.

| Argument | Type | Default | |
|---|---|---|---|
| `:of` | ref | required | a values annotation |
| `:min` | number | none | at least |
| `:max` | number | none | at most |
| `:above` | number | none | greater than |
| `:below` | number | none | less than |

Examples:
- `(values :of cons :above 0.8)`: columns with conservation above 0.8 (case3)

Hints:
- `seqql.ref.kind`: :of takes a values annotation (cons, plddt); for items by value: (pick (items :of A) :where (> (value) 0.8))

### letters
`seq.generator.letters` · generator · result: selection
`(letters L :in Q)`: positions with one of the letters L.

| Argument | Type | Default | |
|---|---|---|---|
| `L` (positional) | text | required | the letters, case-insensitive: C, CHW |
| `:in` | selection | none | where to look; default (all) |

Examples:
- `(and 1brs/A (letters C))`: cysteines of 1brs/A
- `(letters CHW :in p53)`: cysteines, histidines and tryptophans of p53 (case1)

Hints:
- `seqql.arg.missing`: letters needs the letters: (letters C)
- `seqql.limit.positions`: scan less: (and 1brs/A (letters C)) reads 1brs/A only

### positions
`seq.generator.positions` · generator · result: selection · `:where` binds: position
`(positions :in Q :where T)`: positions passing a test.

| Argument | Type | Default | |
|---|---|---|---|
| `:in` | selection | none | where to look; default (all) |
| `:where` | test of a position | none | the test each position passes |

Examples:
- `(positions :in p53 :where (< (index) 10))`: the first ten positions of p53 (case1)

Hints:
- `seqql.limit.positions`: a test runs once per position: give it :in Q, or (and Q (positions :where T))

### search
`seq.generator.search` · generator · result: selection
`(search P :in Q :edits n :classes C :regex b :ends e)`: hits of a pattern, one group each; plain letters (ATG) match themselves, :classes only reads codes such as X or N.

| Argument | Type | Default | |
|---|---|---|---|
| `P` (positional) | text | required | the pattern |
| `:in` | selection | none | where to look; default every unit |
| `:edits` | whole number | `0` | edits allowed |
| `:classes` | text | none | a letter-class table: the pattern's letters are classes (protein, iupac) |
| `:regex` | boolean | `false` | the pattern is a regular expression |
| `:ends` | `local_best_last` \| `local_best_first` \| `all` | `local_best_last` | which ends an approximate search reports |

Examples:
- `(search GXGXXG :in 1brs/A :classes protein)`: a motif
- `(search "[KR]{2}.{10,12}[KR]{3,}" :in p53 :regex true)`: a regular expression, quoted: the bipartite nuclear localization signal (case1)

Hints:
- `seqql.arg.missing`: search needs a pattern: (search GXGXXG)
- `seqql.limit.positions`: search less: (search GXGXXG :in 1brs/A) reads 1brs/A only

### input
`seq.generator.input` · generator · result: selection
`(input name)`: a selection the app passes.

| Argument | Type | Default | |
|---|---|---|---|
| `name` (positional) | text | required | the input's name: selection, focus, structure |

Examples:
- `(input selection)`: the person's selection

Hints:
- `seqql.input.unknown`: inputs are named by the app: (input selection)

### current
`seq.generator.current` · generator · result: selection
`(current)`: inside each, the group being run on.

Examples:
- `(each (items :of regions) (expand (current) 5))`: every region widened by 5 (case1)

Hints:
- `seqql.current`: (current) is only meaningful inside each

## Combine

### and
`seq.combinator.intersect` · combinator · result: selection
`(and Q Q...)`: each group of the first, intersected with the others.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the groups kept |
| `Q…` (positional) | selection | one or more | the selections they are intersected with |

Examples:
- `(and 1brs/A domains)`: positions of 1brs/A inside any domain

Hints:
- `seqql.type`: and mixes a selection and a test: filter with (pick Q :where T)

### or
`seq.combinator.merge` · combinator · result: selection
`(or Q Q...)`: all groups, duplicates dropped.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `Q…` (positional) | selection | one or more | more selections |

Examples:
- `(or (sel "1brs/A[10-20]") (sel "1brs/B[30-35]"))`: two sites, one group each (case7)

### not
`seq.modifier.complement` · modifier · result: selection
`(not Q :in U)`: the positions of U (default all) not in Q.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the positions removed |
| `:in` | selection | none | the positions kept from; default (all) |

Examples:
- `(not domains :in 1brs/A)`: positions of 1brs/A outside every domain

Hints:
- `seqql.arg.count`: not takes one selection; to remove B from A, intersect: (and A (not B))

### range
`seq.modifier.range` · modifier · result: selection
`(range Q part... :numbering N)`: narrow to labels, inclusive.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `part…` (positional) | label | one or more | a label or label range: 10-20, 45A, -5-10 |
| `:numbering` | text | none | the numbering; default the system's |

Examples:
- `(range 1brs/A 10-20 45)`: labels 10 to 20 and 45 of 1brs/A (case7)
- `(range p53 100-120)`: residues 100 to 120 of p53 (case1)

Hints:
- `seqql.arg.order`: positional arguments come before keywords: (range 1brs/A 10-20 :numbering uniprot)
- `seqql.label.unknown`: labels are as the ruler shows them; for indices from 0 use (slice Q start stop)

### slice
`seq.modifier.slice` · modifier · result: selection
`(slice Q start stop...)`: narrow to indices, from 0, end excluded.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `bounds…` (positional) | whole number | one or more | start and stop indices, in pairs: 9 20, or 9 20 30 40 |

Examples:
- `(slice p53 0 50)`: the first 50 positions of p53, as indices (case1)

Hints:
- `seqql.arg.count`: slice takes start-stop pairs: (slice 1brs/A 9 20)
- `seqql.type`: slice takes indices, whole numbers from 0: for labels like 45A use (range 1brs/A 45A)
- `seqql.slice.range`: indices run from 0 to the length; the summary gives each run's slice
- `seqql.slice.clipped`: indices run from 0 to the length; the summary gives each run's slice

### expand
`seq.modifier.expand` · modifier · result: selection
`(expand Q n :left n :right n)`: widen each run.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `n` (positional) | whole number | none | positions added on both sides |
| `:left` | whole number | none | positions added before |
| `:right` | whole number | none | positions added after |

Examples:
- `(and (expand (sel "1brs/A[10-20]") 5) (not (sel "1brs/A[10-20]")))`: around a site, not the site (case7)

Hints:
- `seqql.arg.count`: expand takes one positional after the selection, the width: (expand Q 5)

### onto
`seq.modifier.onto` · modifier · result: selection
`(onto Q R :kind K :keep b)`: the same places in another sequence, codons too; onto an alignment: its columns, :kind row its rows.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the positions carried |
| `R` (positional) | ref | none | a unit, molecule, row or alignment |
| `:kind` | `unit` \| `molecule` \| `row` \| `columns` | none | only for an alignment's rows (row) or columns (columns, the default onto an alignment); without R, every reachable system of the kind; none across a codon mapping |
| `:keep` | boolean | `false` | keep the source runs too |

Examples:
- `(onto (input selection) msa)`: the alignment columns of the selection
- `(onto (sel "protein[5-10]") gene)`: the codons of residues 5-10, no :kind needed (case14)

Hints:
- `seqql.onto.unmapped`: a group with no image is dropped; a row reaches other sequences through its binding
- `seqql.type` (codon): :kind is a kind of sequence (unit, molecule, row, columns), and codon is not one: across a codon mapping leave :kind out, (onto Q gene) gives the codons

### group
`seq.modifier.group` · modifier · result: selection
`(group Q :by G)`: regroup: none (one group), run or system.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `:by` | `none` \| `run` \| `system` | required | none: one group; run: each maximal run; system: one per system |

Examples:
- `(group (letters C :in p53) :by none)`: p53's cysteines as one group (case1)

Hints:
- `seqql.arg.missing`: group needs :by: (group Q :by none)

### partners
`seq.modifier.partners` · modifier · result: selection · `:where` binds: item
`(partners Q :via A :side S :where T)`: the ends of pairs (annotation A) whose other end is in Q.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the partner selection |
| `:via` | ref | required | a pairs annotation |
| `:side` | `any` \| `a` \| `b` | `any` | which ends may be returned: a or b, the pairs' columns (any: both); not a chain, which goes in Q |
| `:where` | test of an item | none | the pairs to use |

Examples:
- `(and 1brs/A (partners 1brs/D :via contacts))`: residues of A in contact with D (case7)

Hints:
- `seqql.ref.kind`: :via takes a pairs annotation (contacts, hbonds)

### each
`seq.modifier.each` · modifier · result: selection
`(each Q E)`: run E per group of Q.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the groups run on |
| `E` (positional) | selection | required | the query run per group, (current) the group |

Examples:
- `(each (items :of regions) (and (current) (letters C)))`: the cysteines of each region, one group each (case1)

### pick
`seq.filter.pick` · filter · result: selection · `:where` binds: group
`(pick Q :where T)`: keep groups passing T.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `:where` | test of a group | required | the test each group passes |

Examples:
- `(pick (items :of domains) :where (matches (label) kinase))`: kinase domains

Hints:
- `seqql.arg.missing`: pick needs a test: (pick Q :where (> (length) 50))

### first
`seq.filter.first` · filter · result: selection
`(first Q n)`: the first n groups.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the selection |
| `n` (positional) | whole number | `1` | how many |

Examples:
- `(first (items :of regions) 2)`: the first two regions (case1)

### ids
`seq.filter.ids` · filter · result: selection
`(ids Q id...)`: keep items with these ids.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | the item groups |
| `id…` (positional) | text | one or more | item ids |

Examples:
- `(ids (items :of domains) PF00069 PF07714)`: two domains by id

Hints:
- `seqql.arg.missing`: ids needs at least one id: (ids (items :of domains) PF00069)

## Tests

### and
`core.logic.and` · core · result: boolean
`(and T T...)`: every test holds; stops at the first that fails.

| Argument | Type | Default | |
|---|---|---|---|
| `T` (positional) | boolean | required | a test |
| `T…` (positional) | boolean | one or more | more tests |

Examples:
- `(items :of regions :where (and (> (length) 30) (matches (label) binding)))`: long binding regions (case1)

Hints:
- `seqql.type`: and mixes a selection and a test: filter with (pick Q :where T)

### or
`core.logic.or` · core · result: boolean
`(or T T...)`: some test holds; stops at the first that does.

| Argument | Type | Default | |
|---|---|---|---|
| `T` (positional) | boolean | required | a test |
| `T…` (positional) | boolean | one or more | more tests |

Examples:
- `(positions :in p53 :where (or (= (letter) C) (= (letter) H)))`: cysteines and histidines of p53 (case1)

Hints:
- `seqql.type`: or mixes a selection and a test: filter with (pick Q :where T)

### not
`core.logic.not` · core · result: boolean
`(not T)`: the test fails; true where its value is missing.

| Argument | Type | Default | |
|---|---|---|---|
| `T` (positional) | boolean | required | a test |

Examples:
- `(items :of regions :where (not (matches (label) binding)))`: regions other than binding ones (case1)

### =
`core.rel.eq` · core · result: boolean
`(= a b)`: equal as JSON values; a missing value equals nothing.

| Argument | Type | Default | |
|---|---|---|---|
| `a` (positional) | value | required | the left side |
| `b` (positional) | value | required | the right side |

Examples:
- `(items :where (= (type) uniprot:Domain))`: every domain item

### !=
`core.rel.neq` · core · result: boolean
`(!= a b)`: not equal as JSON values.

| Argument | Type | Default | |
|---|---|---|---|
| `a` (positional) | value | required | the left side |
| `b` (positional) | value | required | the right side |

Examples:
- `(items :of variants :where (!= (category) benign))`: variants not benign

### <
`core.rel.lt` · core · result: boolean
`(< a b)`: less than; false when a side is missing or not a number.

| Argument | Type | Default | |
|---|---|---|---|
| `a` (positional) | number | required | the left side |
| `b` (positional) | number | required | the right side |

Examples:
- `(positions :in p53 :where (< (index) 10))`: the first ten positions of p53 (case1)

Hints:
- `seqql.type`: < compares numbers: (< (length) 50)

### <=
`core.rel.lte` · core · result: boolean
`(<= a b)`: at most; false when a side is missing or not a number.

| Argument | Type | Default | |
|---|---|---|---|
| `a` (positional) | number | required | the left side |
| `b` (positional) | number | required | the right side |

Examples:
- `(items :of regions :where (<= (length) 40))`: regions of at most 40 positions (case1)

Hints:
- `seqql.type`: <= compares numbers: (<= (length) 50)

### >
`core.rel.gr` · core · result: boolean
`(> a b)`: greater than; false when a side is missing or not a number.

| Argument | Type | Default | |
|---|---|---|---|
| `a` (positional) | number | required | the left side |
| `b` (positional) | number | required | the right side |

Examples:
- `(systems :where (> (length) 100))`: sequences longer than 100 (case1)

Hints:
- `seqql.type`: > compares numbers: (> (length) 50)

### >=
`core.rel.gre` · core · result: boolean
`(>= a b)`: at least; false when a side is missing or not a number.

| Argument | Type | Default | |
|---|---|---|---|
| `a` (positional) | number | required | the left side |
| `b` (positional) | number | required | the right side |

Examples:
- `(pick (letters C) :where (>= (length) 2))`: runs of two or more cysteines (case14)

Hints:
- `seqql.type`: >= compares numbers: (>= (length) 50)

### in
`core.rel.in` · core · result: boolean
`(in x v...)`: x equals one of the v.

| Argument | Type | Default | |
|---|---|---|---|
| `x` (positional) | value | required | the value |
| `v…` (positional) | value | one or more | the values it may equal |

Examples:
- `(positions :in p53 :where (in (letter) C H))`: cysteines and histidines of p53 (case1)

Hints:
- `seqql.arg.missing`: in needs a value and what it may equal: (in (letter) C H)

### matches
`core.str.matches` · core · result: boolean
`(matches s pattern :case b)`: a regular expression, case-insensitive unless :case true.

| Argument | Type | Default | |
|---|---|---|---|
| `s` (positional) | text | required | the text |
| `pattern` (positional) | text | required | an ECMAScript regular expression |
| `:case` | boolean | `false` | case-sensitive |

Examples:
- `(rows :of fam :where (matches (label) HUMAN :case true))`: the human members' rows

Hints:
- `seqql.arg.missing`: matches needs the text and a pattern: (matches (label) kinase)

### index
`seq.position.index` · property · result: whole number · reads: position
`(index)`: its index (0-based, as slice takes it).

Examples:
- `(positions :in 1brs/A :where (< (index) 10))`: the first ten positions of 1brs/A (case7)

Hints:
- `seqql.slot`: index reads a position: use it in (positions :where …)

### letter
`seq.position.letter` · property · result: text · reads: position
`(letter)`: the letter ("" where letters are not held).

Examples:
- `(positions :in p53 :where (= (letter) W))`: tryptophans of p53 (case1)

Hints:
- `seqql.slot`: letter reads a position: use it in (positions :where …)

### label
`seq.position.label` · property · result: text · reads: position
`(label :numbering N)`: its label in the default numbering or :numbering N ("" when none).

| Argument | Type | Default | |
|---|---|---|---|
| `:numbering` | text | none | the numbering; default the system's |

Examples:
- `(positions :in 1brs/A :where (matches (label) "[A-Z]$"))`: positions with an insertion code

Hints:
- `seqql.slot`: label reads a position: use it in (positions :where …)

### value-of
`seq.position.value_of` · property · result: number · reads: position
`(value-of A)`: the value of values annotation A here; NaN when none.

| Argument | Type | Default | |
|---|---|---|---|
| `A` (positional) | ref | required | a values annotation |

Examples:
- `(positions :in model :where (> (value-of plddt) 70))`: confident positions (case8)

Hints:
- `seqql.slot`: value-of reads a position: use it in (positions :where …)
- `seqql.ref.kind`: value-of takes a values annotation (cons, plddt)

### ref
`seq.system.ref` · property · result: text · reads: system
`(ref)`: the system's ref.

Examples:
- `(systems :where (matches (ref) ^1brs/))`: the units of 1brs (case7)

Hints:
- `seqql.slot`: ref reads a system: use it in (systems :where …) or (rows :where …)

### kind
`seq.system.kind` · property · result: text · reads: system
`(kind)`: the system's kind: unit, molecule, row or columns.

Examples:
- `(pick (all) :where (= (kind) row))`: every alignment row (case3)

Hints:
- `seqql.slot`: kind reads a system: use it in (systems :where …) or (rows :where …)

### label
`seq.system.label` · property · result: text · reads: system
`(label)`: the system's label.

Examples:
- `(rows :of fam :where (matches (label) HUMAN))`: the human members' rows

Hints:
- `seqql.slot`: label reads a system: use it in (systems :where …) or (rows :where …)

### length
`seq.system.length` · property · result: whole number · reads: system
`(length)`: positions in the system.

Examples:
- `(systems :where (> (length) 100))`: sequences longer than 100 (case1)

Hints:
- `seqql.slot`: length reads a system: use it in (systems :where …) or (rows :where …)

### origin
`seq.system.origin` · property · result: text · reads: system
`(origin)`: the system's origin.

Examples:
- `(systems :where (= (origin) sequence))`: units that are their molecule's own sequence (case1)

Hints:
- `seqql.slot`: origin reads a system: use it in (systems :where …) or (rows :where …)

### score
`seq.system.score` · property · result: number · reads: system
`(score)`: a unit's score.

Examples:
- `(systems :where (>= (score) 0.5))`: units scoring at least 0.5

Hints:
- `seqql.slot`: score reads a system: use it in (systems :where …) or (rows :where …)

### alphabet
`seq.system.alphabet` · property · result: text · reads: system
`(alphabet)`: the system's alphabet: protein, dna, rna or generic.

Examples:
- `(systems :where (= (alphabet) dna))`: DNA sequences (case14)

Hints:
- `seqql.slot`: alphabet reads a system: use it in (systems :where …) or (rows :where …)

### field
`seq.system.field` · property · result: value · reads: system
`(field name :of K :ns N)`: another field: :of a row table, :ns app fields.

| Argument | Type | Default | |
|---|---|---|---|
| `name` (positional) | text | required | the field: molecule, complex, alignment, row, or a column of :of or :ns |
| `:of` | ref | none | a rows key table |
| `:ns` | text | none | the namespace of an app's custom fields |

Examples:
- `(systems :where (= (field auth_asym_id :ns mvs) A))`: units of chain A, by the app's field (case7)

Hints:
- `seqql.slot`: field reads a system: use it in (systems :where …) or (rows :where …)
- `seqql.field.unknown`: a system's own fields: molecule, complex, alignment, row; others need :of or :ns

### id
`seq.item.id` · property · result: text · reads: item
`(id)`: the item's id.

Examples:
- `(pick (items :of domains) :where (in (id) PF00069 PF07714))`: two domains by id

Hints:
- `seqql.slot`: id reads an item: use it in (items :where …) or (partners :where …)

### label
`seq.item.label` · property · result: text · reads: item
`(label)`: the item's label.

Examples:
- `(items :of regions :where (matches (label) binding))`: p53's DNA-binding region (case1)

Hints:
- `seqql.slot`: label reads an item: use it in (items :where …) or (partners :where …)

### category
`seq.item.category` · property · result: text · reads: item
`(category)`: the item's category.

Examples:
- `(partners 1brs/D :via contacts :where (= (category) hbond))`: ends hydrogen-bonded to 1brs/D

Hints:
- `seqql.slot`: category reads an item: use it in (items :where …) or (partners :where …)

### value
`seq.item.value` · property · result: value · reads: item
`(value)`: the item's value.

Examples:
- `(pick (items :of domains) :where (> (value) 0.8))`: domains valued above 0.8

Hints:
- `seqql.slot`: value reads an item: use it in (items :where …) or (partners :where …)

### type
`seq.item.type` · property · result: text · reads: item
`(type)`: its annotation's semantic type.

Examples:
- `(items :where (= (type) uniprot:Domain))`: every domain item

Hints:
- `seqql.slot`: type reads an item: use it in (items :where …) or (partners :where …)

### annotation
`seq.item.annotation` · property · result: text · reads: item
`(annotation)`: its annotation's ref.

Examples:
- `(items :where (= (annotation) variants))`: the items of variants

Hints:
- `seqql.slot`: annotation reads an item: use it in (items :where …) or (partners :where …)

### length
`seq.item.length` · property · result: whole number · reads: item
`(length)`: positions in the item.

Examples:
- `(items :of regions :where (> (length) 50))`: regions longer than 50 (case1)

Hints:
- `seqql.slot`: length reads an item: use it in (items :where …) or (partners :where …)

### field
`seq.item.field` · property · result: value · reads: item
`(field name)`: any other column.

| Argument | Type | Default | |
|---|---|---|---|
| `name` (positional) | text | required | the column or field: tooltip, fragments, evalue |

Examples:
- `(items :of hits :where (< (field evalue) 0.001))`: hits with an e-value below 0.001

Hints:
- `seqql.slot`: field reads an item: use it in (items :where …) or (partners :where …)
- `seqql.field.unknown`: an annotation's own columns are listed when the name is unknown

### length
`seq.group.length` · property · result: whole number · reads: group
`(length)`: positions in the group.

Examples:
- `(pick (letters C) :where (>= (length) 2))`: runs of two or more cysteines (case14)

Hints:
- `seqql.slot`: length reads a group: use it in (pick Q :where …)

### runs
`seq.group.runs` · property · result: whole number · reads: group
`(runs)`: the group's run count.

Examples:
- `(pick (items :of domains) :where (> (runs) 1))`: domains in more than one piece

Hints:
- `seqql.slot`: runs reads a group: use it in (pick Q :where …)

### score
`seq.group.score` · property · result: number · reads: group
`(score)`: a hit's edit count; NaN otherwise.

Examples:
- `(pick (search GXGXXG :edits 1) :where (= (score) 0))`: exact hits only

Hints:
- `seqql.slot`: score reads a group: use it in (pick Q :where …)

### overlaps
`seq.group.overlaps` · property · result: boolean · reads: group
`(overlaps Q)`: shares a position with Q.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | computed once, not per group |

Examples:
- `(pick sites :where (overlaps (pick domains :where (matches (label) SH2))))`: sites in an SH2 domain (a bare name selects all of the annotation)
- `(pick (items :of variants) :where (overlaps (sel "domains[#PF00069]")))`: variants in a domain, by its id

Hints:
- `seqql.slot`: overlaps reads a group: use it in (pick Q :where …)

### within
`seq.group.within` · property · result: boolean · reads: group
`(within Q n)`: within n positions of Q, in the same system.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | computed once, not per group |
| `n` (positional) | whole number | required | the distance in positions |

Examples:
- `(pick (items :of sites) :where (within (sel "1brs/A[40-90]") 5))`: sites near a stretch

Hints:
- `seqql.slot`: within reads a group: use it in (pick Q :where …)

### count
`seq.group.count` · property · result: whole number · reads: group
`(count Q)`: how many of Q's positions it holds.

| Argument | Type | Default | |
|---|---|---|---|
| `Q` (positional) | selection | required | computed once, not per group |

Examples:
- `(pick (items :of domains) :where (>= (count (letters C)) 4))`: domains with four or more cysteines

Hints:
- `seqql.slot`: count reads a group: use it in (pick Q :where …)

## Issue codes

| Code | Severity |
|---|---|
| `seqql.read.unclosed` | error |
| `seqql.read.extra_close` | error |
| `seqql.read.trailing` | error |
| `seqql.read.string` | error |
| `seqql.read.reserved` | error |
| `seqql.read.keyword` | error |
| `seqql.read.empty` | error |
| `seqql.read.depth` | error |
| `seqql.read.size` | error |
| `seqql.read.head` | error |
| `seqql.reply.empty` | error |
| `seqql.reply.nested` | error |
| `seqql.reply.none` | error |
| `seqql.reply.unclosed` | error |
| `seqql.reply.unopened` | error |
| `seqql.symbol.unknown` | error |
| `seqql.symbol.duplicate` | error |
| `seqql.symbol.namespace` | error |
| `seqql.symbol.unavailable` | error |
| `seqql.alias.collision` | error |
| `seqql.alias.reserved` | error |
| `seqql.arg.count` | error |
| `seqql.arg.missing` | error |
| `seqql.arg.order` | error |
| `seqql.arg.repeated` | error |
| `seqql.arg.unknown` | error |
| `seqql.type` | error |
| `seqql.slot` | error |
| `seqql.sel.syntax` | error |
| `seqql.sel.mixed` | error |
| `seqql.macro.depth` | error |
| `seqql.current` | error |
| `seqql.field.unknown` | error |
| `seqql.slice.range` | error |
| `seqql.ref.unknown` | error |
| `seqql.ref.kind` | error |
| `seqql.numbering.unknown` | error |
| `seqql.numbering.ambiguous` | error |
| `seqql.classes.unknown` | error |
| `seqql.input.unknown` | error |
| `seqql.limit.positions` | error |
| `seqql.limit.runs` | error |
| `seqql.limit.groups` | error |
| `seqql.limit.time` | error |
| `seqql.aborted` | error |
| `seqql.lazy` | warning |
| `seqql.letters.columns` | warning |
| `seqql.label.unknown` | warning |
| `seqql.ids.unknown` | warning |
| `seqql.onto.kind` | warning |
| `seqql.onto.ambiguous` | warning |
| `seqql.onto.unmapped` | warning |
| `seqql.field.near` | warning |
| `seqql.slice.clipped` | warning |
| `seqql.empty` | info |
