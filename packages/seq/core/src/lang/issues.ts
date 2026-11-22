/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `seq.lang.*` codes: what the description cannot express (_spec/README.md, Issues). Each with its severity, what
// it means and the usual next step, as the tree core's catalogs are; a test fails on a code a source uses and the
// catalog lacks.

import { defineIssues } from '@molstar/core/repr/tree';

export const LangIssues = defineIssues('seq.lang', {
  'seq.lang.version.unsupported': {
    severity: 'error',
    description: 'A spec is written in a newer major version of the language than this one reads',
    hint: 'Read it with a newer version of the library, or migrate the spec down',
  },
  'seq.lang.id.required': {
    severity: 'error',
    description: 'A node of the spec has no id (an extension kind’s too), so it cannot be addressed or patched',
    hint: 'Give the node an id unique among its siblings',
  },
  'seq.lang.container.id': {
    severity: 'error',
    description:
      'A container (`data`, `views`, `refs`, `snapshots`, `navigators`, `compounds`) has an id other than its kind',
    hint: 'The containers have reserved ids',
  },
  'seq.lang.source.index': {
    severity: 'error',
    description: 'A source with an inline `value` has an `index`, which only a file read by region has',
    hint: 'Drop the `index`, or give a `uri` instead of the `value`',
  },
  'seq.lang.source.content': {
    severity: 'error',
    description:
      'A source has both a `uri` and an inline `value`, or neither (nor a host’s `ref`), or a `ref` beside a `value`',
    hint: 'Give one of `uri` and `value`; a host’s `ref` beside a `uri`, or alone; inline data is never a host’s',
  },
  'seq.lang.snapshot.content': {
    severity: 'error',
    description: 'A story step has both a `patch` and a `state`, or neither',
    hint: 'Give one: a patch over the spec (or the step it starts `from`), or a full state',
  },
  'seq.lang.sequence.content': {
    severity: 'error',
    description: 'A sequence has both `letters` and a `source`, or neither (an override under an assembly: either)',
    hint: 'Give one of `letters` and `source`; an override reads its assembly’s',
  },
  'seq.lang.alignment.content': {
    severity: 'error',
    description: 'An alignment has both inline `rows` and a `source`, or neither',
    hint: 'Give one of `rows` and `source`',
  },
  'seq.lang.alignment.width': {
    severity: 'error',
    description: 'An inline alignment’s rows are not all as wide, or its labels are not one a row',
    hint: 'Pad every row with gaps to the alignment’s width; give a label a row',
  },
  'seq.lang.sequence.letters': {
    severity: 'error',
    description: 'Inline letters of a sequence are not letters of its alphabet',
    hint: 'Protein: letters, `*` and `-`; dna: ACGTN, the ambiguity codes and `-`; rna: the same with U',
  },
  'seq.lang.ranges.content': {
    severity: 'error',
    description: 'A ranges node has both `rows` and a `source`, or neither',
    hint: 'Give one of `rows` and `source`',
  },
  'seq.lang.ranges.interval': {
    severity: 'error',
    description: 'A row of ranges starts before position 0, or ends before it starts on a sequence that is not cyclic',
    hint: 'Ranges are 0-based and half-open: 0 <= start <= end; only a cyclic sequence has ranges that wrap',
  },
  'seq.lang.ranges.strand': {
    severity: 'error',
    description: 'A strand of a row of ranges is neither `+` nor `-`',
    hint: 'The strand column holds `+` or `-`',
  },
  'seq.lang.ranges.bounds': {
    severity: 'error',
    description: 'A row of ranges ends past the inline letters of its target sequence',
    hint: 'A range ends at most at the length of the sequence',
  },
  'seq.lang.pairs.interval': {
    severity: 'error',
    description: 'An end of a pair starts before position 0, or does not end after it starts',
    hint: 'Ends are 0-based and half-open, at least one position, and never wrap: 0 <= start < end',
  },
  'seq.lang.pairs.bounds': {
    severity: 'error',
    description: 'An end of a pair ends past the inline letters of its sequence',
    hint: 'An end ends at most at the length of its sequence',
  },
  'seq.lang.grid.length': {
    severity: 'error',
    description: 'A grid’s values are not one per pair of positions of its two sequences',
    hint: 'Give |x| × |y| values, x-major; NaN is missing',
  },
  'seq.lang.values.content': {
    severity: 'error',
    description: 'A values node has both `values` and a `source`, or neither, or a `field` without a `source`',
    hint: 'Give one of `values` and `source`; `field` names the source’s column of values',
  },
  'seq.lang.values.length': {
    severity: 'error',
    description: 'Inline values are not one per position of the inline letters of their target sequence',
    hint: 'Give as many values as the sequence has letters; NaN is missing',
  },
  'seq.lang.numbering.content': {
    severity: 'error',
    description:
      'A numbering’s params do not fit its method: `labels` with both `labels` and a `source`, or neither, or a `field` without a `source`; `offset` without `by`; `segments` without `segments`; or a param another method takes',
    hint: '`labels` takes one of `labels` and `source` (`field` names the source’s column); `offset` takes `by`; `segments` takes `segments`',
  },
  'seq.lang.numbering.segments': {
    severity: 'error',
    description:
      'A numbering’s run ends before it starts or overlaps another, or one with an insertion code holds more than one position',
    hint: 'Give each run `[start, end)` apart from the others; an insertion code (`52A`) numbers one position',
  },
  'seq.lang.numbering.id': {
    severity: 'error',
    description: 'A numbering’s id is a built-in numbering’s (`index`, `ordinal`): forms would never read it',
    hint: 'Name it by its scheme (`kabat`, `imgt`, `auth`)',
  },
  'seq.lang.numbering.length': {
    severity: 'error',
    description: 'Inline labels are not one a position of the inline letters of their sequence',
    hint: 'Give as many labels as the sequence has letters; `null` is no label there',
  },
  'seq.lang.numbering.duplicate': {
    severity: 'warning',
    description: 'A label is given at two positions: it names the first',
    hint: 'Give each position its own label (an insertion code: `100A`), or `null`',
  },
  'seq.lang.map.params': {
    severity: 'error',
    description: 'A map has params its method does not take, or lacks one it needs',
    hint: '`offset` takes `by`; `codon` takes `exons` and `strand`; `segments` takes `segments`; `structure` takes `segments` (with `provenance`) and `status`',
  },
  'seq.lang.map.segments': {
    severity: 'error',
    description:
      'A map’s run holds no position, overlaps another on `from` (a structure map’s on `to` too), or has a strand other than `+` and `-`; or a status row ends before it starts or says other than `conflict` and `unobserved`',
    hint: 'Give each run a length of at least 1, apart from the others on `from`; a status row says `conflict` or `unobserved`',
  },
  'seq.lang.map.alphabet': {
    severity: 'error',
    description: 'A codon map does not go from a protein to a dna or rna sequence',
    hint: 'A codon map is a protein on its gene',
  },
  'seq.lang.map.exons': {
    severity: 'error',
    description: 'An exon of a codon map ends before it starts, or overlaps or precedes the exon before it',
    hint: 'List the exons in genomic order, each from its start to its end, none overlapping',
  },
  'seq.lang.ref.unresolved': {
    severity: 'warning',
    description:
      'A user-defined ref names nothing: the leading ids of its address name no node of the spec (the rest of the spec opens)',
    hint: 'Point the ref at an address of the spec, or remove it',
  },
  'seq.lang.ref.cycle': {
    severity: 'warning',
    description: 'A user-defined ref leads back to itself through other refs, or through more refs than are followed',
    hint: 'Point one of the refs at what it names, not at another ref',
  },
  'seq.lang.ref.target': {
    severity: 'error',
    description: 'A reference names no node of the spec, or one of a kind the param does not take',
    hint: 'Name a node by its address, like /data/p53',
  },
  'seq.lang.layer.data': {
    severity: 'error',
    description:
      'A layer’s data does not fit its type: `ranges` draws ranges, `values` draws values, `letters` draws letters',
    hint: 'Give data of the kind the type draws',
  },
  'seq.lang.layer.fallback': {
    severity: 'warning',
    description: 'A layer of an extension’s type has no `fallback` (read as `none`)',
    hint: 'Give a built-in type, or `none`',
  },
  'seq.lang.section.range': {
    severity: 'error',
    description:
      'The range of a section ends before it starts on a sequence that is not cyclic, or runs past the inline letters of its system',
    hint: 'A range is [start, end) within the system',
  },
  'seq.lang.theme.form': {
    severity: 'error',
    description:
      'A channel slot holds a node that is no form (a palette, a clause or a measure where a form goes), or a palette or clause slot holds a form',
    hint: 'A slot holds a string, a number, or a `uniform`, `field`, `rules`, `scheme` or `letters` form',
  },
  'seq.lang.theme.color': {
    severity: 'error',
    description: 'A color in a theme is no color: not `#rgb`/`#rrggbb`(aa), `token:<name>` or `light-dark(a, b)`',
    hint: 'Write `#d62728`, `token:range` or `light-dark(#333333, #dddddd)`',
  },
  'seq.lang.theme.channel': {
    severity: 'error',
    description:
      'A value not of the channel’s kind: a color in `size`, a number in `color`, a word the channel does not take',
    hint: 'Color channels take colors, `size` and `opacity` numbers, `line`, `arrow` and `shape` their words',
  },
  'seq.lang.theme.type': {
    severity: 'warning',
    description:
      'A clause compares across types (text with a number, or `<` with text): it never holds. Text compares with text, numbers with numbers',
    hint: 'Compare a number column with numbers, a text column with text',
  },
  'seq.lang.theme.stops': {
    severity: 'error',
    description: 'A palette’s edges or stops, or a domain’s ends, are out of order',
    hint: 'Edges and stops ascend; a domain is `[low, high]`',
  },
  'seq.lang.theme.scale': {
    severity: 'error',
    description: 'A theme names a scale by an address that names no `scale` node',
    hint: 'Name a scale under `views`, like /views/conservation',
  },
  'seq.lang.theme.depth': {
    severity: 'error',
    description: 'Theme forms nest past the bound (`SpecBounds.themeDepth`, 4)',
    hint: 'Flatten the rules: a rule’s value a value, or one form',
  },
  'seq.lang.section.follow': {
    severity: 'error',
    description: 'A section follows itself, directly or through other sections',
    hint: 'Following is a chain that ends at a section that follows none',
  },
  'seq.lang.tooltip.content': {
    severity: 'error',
    description:
      'A tooltip names `on` beside a table, or neither: it is by data (`on`) or by range (`target` with `rows` or a `source`)',
    hint: 'Give `on` alone, or a `target` with one of `rows` and `source`',
  },
  'seq.lang.tooltip.on': {
    severity: 'warning',
    description:
      'A tooltip’s `on`, or an address of its table, names nothing the view draws: that row is inert, the rest shows',
    hint: 'Name a data node, a layer, or an address that resolves in the spec',
  },
  'seq.lang.tooltip.field': {
    severity: 'warning',
    description: 'A tooltip’s template names a field its table lacks: no row of it gives a text',
    hint: 'Name a column of the table (`{}` reads `field`, `tooltip` when absent)',
  },
  'seq.lang.tooltip.schema': {
    severity: 'error',
    description: 'A tooltip table’s columns do not fit its `schema` (none of the columns the schema locates by)',
    hint: 'Say the schema the table is written in (`ranges`: `start` and `end`; `address`: `address`; MVS’s by name)',
  },
  'seq.lang.tooltip.category': {
    severity: 'info',
    description:
      'A tooltip reads a CIF category, which only a structure’s source holds and this reader does not read yet: it shows nothing',
    hint: 'Give the table inline or from a source of its own until the structure runtime reads categories',
  },
  'seq.lang.tooltip.text': {
    severity: 'error',
    description:
      'A tooltip’s template is not an f-string this reader formats, or is past the bound (`SpecBounds.tooltipText`, 1,000)',
    hint: 'Write fields as `{name}` or `{}` (`{{` for a brace), with a format a value takes; keep it short',
  },
  'seq.lang.navigator.bands': {
    severity: 'error',
    description: 'A navigator draws `bands` and names no annotation of bands',
    hint: 'Name the bands (`bands`: a ranges node with a `stain` field), or draw another `content`',
  },
  'seq.lang.source.host_only': {
    severity: 'info',
    description:
      'A source is named by a host’s `ref` alone: without that host it is unavailable, so the spec is not portable (EC-185)',
    hint: 'Give the data’s `uri` beside the `ref` where there is one',
  },
  'seq.lang.modified.order': {
    severity: 'error',
    description: 'A sequence’s modified positions are not in order (rows at one position are alternatives)',
    hint: 'List the rows by position, ascending',
  },
  'seq.lang.modified.bounds': {
    severity: 'error',
    description: 'A modified position starts before 0, or it (with its `length`) ends past the inline letters',
    hint: 'A row’s position and length lie within the sequence',
  },
  'seq.lang.modified.length': {
    severity: 'error',
    description: 'A modified row’s `length` is less than 1, or its component runs into the next row’s position',
    hint: 'A component covers its `length` positions (1 when absent), up to the next row',
  },
  'seq.lang.molecule.content': {
    severity: 'error',
    description: 'A molecule has both `residues` and a `source`, or neither, or `atoms` without inline `residues`',
    hint: 'Give one of `residues` and `source`; its atoms beside its residues',
  },
  'seq.lang.molecule.atoms': {
    severity: 'error',
    description: 'A molecule’s atoms name a residue it does not have, or their residues decrease',
    hint: 'Each atom’s `residue` is a position of `residues`, the column non-decreasing: a residue’s atoms are a run',
  },
  'seq.lang.molecule.hydrogen': {
    severity: 'error',
    description: 'A molecule’s atoms hold a hydrogen: they are heavy atoms only (EC-64)',
    hint: 'Leave the hydrogens out: they are often absent or model-dependent, so lengths would not line up',
  },
  'seq.lang.molecule.identifier': {
    severity: 'error',
    description: 'An identifier’s namespace is not an open lowercase word (`ccd`, `smiles`, a host’s own)',
    hint: 'Write the namespace lowercase: a letter, then letters, digits, `_` or `-`',
  },
  'seq.lang.molecule.form': {
    severity: 'error',
    description: 'A block is drawn as `atoms`, and its system has no atoms (a sequence, or a molecule without them)',
    hint: 'Give the molecule its `atoms`, or draw the block `whole` or by `residues`',
  },
  'seq.lang.unit.range': {
    severity: 'error',
    description: 'A unit’s range ends before it starts (a unit never wraps, EC-13), or runs past its molecule',
    hint: 'A unit’s range is [start, end) within the positions of what it instantiates',
  },
  'seq.lang.unit.observed': {
    severity: 'error',
    description: 'A unit’s observed runs are not intervals in order on its own positions',
    hint: 'List the modeled runs in order on the unit’s own positions (from 0), each [start, end), none overlapping',
  },
  'seq.lang.link.target': {
    severity: 'error',
    description:
      'A pairs node under `data` names no `a`, or links (a child of a complex or a molecule) name an `a` or `b`',
    hint: 'Under `data` give `a`; under a complex or a molecule, the parent is both ends',
  },
  'seq.lang.link.unit': {
    severity: 'error',
    description:
      'A link of a complex names no unit of it at an end, or a pairs node not under a complex has a unit column',
    hint: 'Give each row its `a_unit` and `b_unit`, ids of the complex’s units; elsewhere the ends are the node’s',
  },
  'seq.lang.link.bounds': {
    severity: 'error',
    description: 'A link’s end lies outside its unit’s range, or past its molecule’s residues',
    hint: 'An end is a position of its unit (or of the molecule it is under)',
  },
  'seq.lang.link.position': {
    severity: 'error',
    description: 'A link’s end (under a complex or a molecule) is more than one position',
    hint: 'A link joins one position to one position; a site or a contact of several is a pairs node under `data`',
  },
  'seq.lang.assembly.content': {
    severity: 'error',
    description: 'An assembly has no `sequences` or `source`, or `names` alone',
    hint: 'Give one',
  },
  'seq.lang.sequence.alphabet': {
    severity: 'error',
    description: 'A sequence has no `alphabet`',
    hint: 'Give one',
  },
  'seq.lang.block.scale': {
    severity: 'error',
    description: 'A block’s `scale` is 0',
    hint: 'Give a positive number',
  },
  'seq.lang.bridge.sides': {
    severity: 'error',
    description: 'A bridge joins a section to itself',
    hint: 'Name two sections',
  },
  'seq.lang.section.include': {
    severity: 'warning',
    description: '`include` on a section of no assembly is ignored',
    hint: 'Drop it',
  },
  'seq.lang.include.pattern': {
    severity: 'error',
    description: 'An `include` pattern outside the language',
    hint: '`^`, `.`, `[…]`, `?`, `*`, `+`, `$`',
  },
  'seq.lang.pairs.content': {
    severity: 'error',
    description: 'A pairs node has not one of `rows` and `source`',
    hint: 'Give one',
  },
  'seq.lang.grid.content': {
    severity: 'error',
    description: 'A grid has not one of `values` and `source`',
    hint: 'Give one',
  },
  'seq.lang.matrix.content': {
    severity: 'error',
    description: 'A matrix has not one of its cells and `source`',
    hint: 'Give one',
  },
  'seq.lang.data.names': {
    severity: 'error',
    description: 'Inline rows over an assembly name no sequence, or rows over a sequence name one',
    hint: 'Name them over an assembly only, one a row',
  },
  'seq.lang.link.occupancy': {
    severity: 'error',
    description: 'A link’s occupancy is not between 0 and 1',
    hint: 'Occupancy is the share of copies a link is in: 0 to 1 (NaN: 1)',
  },
});

/** The bounds' codes (`bounds.ts`): what a spec sizes, refused before it is read (_spec/README.md, Language, Bounds). */
export const BoundIssues = defineIssues('seq.lang.bound', {
  'seq.lang.bound.size': {
    severity: 'error',
    description:
      'A spec is larger decoded than the configuration allows (`bounds.maxBytes`, 64 MB); it is not read further',
    hint: 'Put the data in sources the spec names, or raise the bound where the application composes it',
  },
  'seq.lang.bound.depth': {
    severity: 'error',
    description:
      'The spec nests deeper than the configuration allows (`bounds.maxDepth`, 256 levels of nodes and values)',
    hint: 'Flatten the tree',
  },
  'seq.lang.bound.shape': {
    severity: 'error',
    description: 'A node built in code is not a node: no record, no kind, or a children list that is no list',
    hint: 'A node is a record with a `kind`; `children` is a list of nodes',
  },
  'seq.lang.bound.query_depth': {
    severity: 'error',
    description: 'A query form nests deeper than the configuration allows (`bounds.maxQueryDepth`, 64)',
    hint: 'Name the nodes with refs instead of nested forms',
  },
  'seq.lang.bound.column': {
    severity: 'error',
    description: 'A column is longer than the configuration allows (`bounds.maxColumnLength`, 2^31 − 1)',
    hint: 'Give the rows in a source, read by region',
  },
  'seq.lang.bound.issues': {
    severity: 'warning',
    description:
      'A code was reported more often than the configuration keeps (`bounds.maxPlaces` places a code, `bounds.maxNodeIssues` a node); the rest are counted, not listed',
    hint: 'Fix the first places; the rest are likely the same cause',
  },
});
