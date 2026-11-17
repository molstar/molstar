# Addressing: everything in a view has an address

## Intent

The author, 2026-11-17: "markdown links etc. . this is about making anything in the view addressable .. so we should
have that ... we should [have] storybooks that showcase the addressing stuff"; then "we should be able to address
"rects", and the syntax for the address should be path/(query :key) for make this generic. meso uses this schema too.
and for identifers: use hierarchical addressig? ipv6 llike? and we should have "user-defined refs"".

A paper's text, a notebook, a chat with an agent, a bookmark, another viewer: each needs to point at something (a node,
a residue, a range in a numbering, a feature, a row, a rectangle of a dot plot or an alignment, a track, a snapshot) and
have the view act on it. One grammar serves all of them, the same as meso\*'s.

## Precedents

- meso\*: local ids, an address the path of ids from the root (as an IPv6 address is a path of groups), `/x/*` for a
  subtree, query forms as segments (`/cytoplasm/(population :ingredient catalase)`), one word wherever an address
  appears, bounded selection, `ref` and `def` params and definition containers ([sources/meso](../sources/meso/README.md);
  its core spec's "Ids and addresses", "Query segments", "References", "Definitions").
- The previous implementation: global refs and derived refs (`1brs/A`, `msa/3`); markdown links in Mol\*'s grammar
  (`!key=value&…`) carrying Mol\*'s and Seq\*'s keys ([sources/seqstar](../sources/seqstar/README.md)).

## Decision

### Ids: local and hierarchical

- **An id is local**: unique among its siblings, not across the spec. A node's address is the path of ids from the root:
  `/data/p53/features`. An id is one segment: a letter, digit or `_`, then letters, digits, `_`, `.`, `+`, `-`. Why: short,
  composable specs (two `features` under two sequences), and an address is a bare word with nothing to quote.
- **Prefixes name subtrees**: an address is a prefix of every address below it; `/data/p53/*` names the subtree.
- **No elision**: an address of ids is exact; searching at any depth is the query language's job, not an address's.
- **Nodes made from data have ids from their keys**: a structure's chains (`/data/1tup/A`), an alignment's rows
  (`/data/globins/rows/HBA_HUMAN`), an annotation's items where they have ids. This replaces derived refs.
- **Names that are not ids** (an HLA contig `HLA-A*01:01:01:01`, EC-200) are reached by a form,
  `(sequence :name "HLA-A*01:01:01:01")`, never by loosening the id pattern.

### Addresses: `path/(query :key value)`

- **The grammar is meso\*'s**:

  ```
  address = "/" | ("/" step)+ ["/*"]
  step    = id | form
  form    = "(" kind (":" key value)* ")"
  value   = atom | form
  ```

  An address of ids names one node; with forms it names a set (each result an exact address). In a spec, a link, a
  command or agent text it is one word, read and printed canonically (`parse(print(a)) = a`).
- **Forms are the query language's**, so an address and a query are one language: a form selects among what the
  prefix names. Over nodes, by kind and params (`/views/main/(section :of /data/p53)`); below a node that has
  positions, items or rows, by the forms of what it holds:

  | Address | Names |
  |---|---|
  | `/data/p53/(range :from 120 :to 130 :numbering uniprot)` | residues 120 to 130 in UniProt's numbering |
  | `/data/1tup/A/(at :label 94A :numbering auth)` | one residue by its author number and insertion code |
  | `/data/p53/features/(item :id DNA-binding)` | one feature |
  | `/data/p53/features/(item :type domain)` | every domain |
  | `/data/globins/rows/(row :species human)` | rows by a field |
  | `/data/globins/(rect :rows (rows :from HBA_HUMAN :to HBB_HUMAN) :columns (range :from 30 :to 60))` | a block of an alignment |
  | `/data/trna/pairs/(rect :a (range :from 10 :to 25) :b (range :from 40 :to 55))` | a rectangle of a dot plot or grid |
  | `/data/meth/(rect :states (states :in [s1 s2]) :positions (range :from 7610000 :to 7620000))` | a block of a matrix |
  | `/views/main/(viewport :focus /data/p53/(range :from 90 :to 300))` | a view opened on a focus |
  | `/snapshots/step3` | a snapshot |

- **Rects** are first-class: a rectangle is two ranges over two axes (positions × positions for a grid or dot plot,
  rows × columns for an alignment, states × positions for a matrix), selected, highlighted, focused and linked like a
  range; a dot-plot drag selects one (EC-145).
- **Bounded**: selection takes a limit and a depth, as meso\*'s does; more matches is an issue, never a truncation.
- **Positions in labels**: forms take labels of a named numbering, resolved at the edge (EC-1 to EC-8); an address that
  no longer resolves is reported with what it named, never silently re-bound.

### User-defined refs

- **Authors name things**: a `refs` container at the root holds named refs, each an id and an address
  (`{ kind: "ref", id: "dbd", params: { to: "/data/p53/(range :from 94 :to 292 :numbering uniprot)", color, label } }`).
  A ref's address is `/refs/dbd`, usable wherever an address is: in layers, focus, commands, links and other refs.
- **Bookmarks and location sets are refs**: adding a bookmark adds a ref (one patch, one undo step, EC-162); a set of
  them is a ref to several addresses.
- **Refs follow edits**: a move or a rename rewrites what names it; a ref whose target no longer resolves is reported,
  and the rest of the spec opens.
- **A `ref` param type** (meso\*'s): a param that names a node or location holds an address word, checked in one
  resolve pass over one index of the spec.

### Actions and links

- **Commands take addresses**: focus, reveal, select, highlight (a named channel), hover, step a snapshot. A link, a
  menu, a key, an agent and a test drive the view the same way.
- **Links** carry an action and an address in Mol\*'s link grammar (`!seq-focus=/refs/dbd`), so one href can carry
  Mol\*'s and Seq\*'s keys, each viewer running the ones it knows; a link's hover uses its own channel (EC-182).
- **Copy link**: the view gives the address of whatever is under the pointer or selected, as the shortest stable form
  (a ref if one names it, else ids and forms).

## Tensions

- **Local ids vs stable identity**: an address changes when a node moves; a move is an operation that rewrites what
  names it (meso\*'s resolution). Accepted for short, composable specs.
- **The id pattern vs real names** (EC-200): kept strict; such names go through a `:name` form.

## Consequences

- Address parsing, printing and selection are conformance-tested with the numbering cases (insertion codes, IMGT
  order, cyclic ranges, reversed pieces) and with rects on every geometry that has two axes.
- Agent tools (deferred) read and write through addresses and expressions; the addressing grammar is their vocabulary.
- The Storybook showcases addressing: a text with links that focus, highlight and select in the view and in the
  structure; copy link; one location in several numberings; a rect from a dot plot and from an alignment; refs as
  bookmarks.
