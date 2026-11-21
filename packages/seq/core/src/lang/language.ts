/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqViewSpec, described as data on `@molstar/core/repr` (_spec/README.md, Language): the first slice's kinds. What
// the description cannot say (a sequence has letters or a source, a codon map links a protein to a gene) is the checks
// of `check.ts`, in the `seq.lang.*` codes.

import { patchLanguage } from '@molstar/core/repr/patch';
import { defineLanguage, type Language, type ParamType, type ValidateOptions } from '@molstar/core/repr/schema';
import {
  BlockForm,
  LettersModified,
  LinkColumns,
  ModifiedParam,
  MoleculeKinds,
  SectionClasses,
  SystemKinds,
} from './molecule-language.js';
import { ScaleParamsDescription, ThemeLanguage, ThemeLanguageId } from './theme-language.js';

/**
 * The version of the language this description is (the package's `LanguageVersion` is the same). Pre-release (`0.x`):
 * a change moves it with no migration (_spec/README.md, Versions).
 */
export const SeqLanguageVersion = '0.5.0';

/**
 * What `pairs.a`/`b`, `ranges.target` and `map.to` take: a system, or an alignment by its columns system (its address;
 * _spec/README.md, Language, Systems).
 */
const ColumnsTargetKinds = [...SystemKinds, 'alignment'] as const;

/**
 * An extension's name, `<id>:<name>` (_spec/README.md, Language, Extensions): its id the namespace, never `seq`
 * (reserved: the check says so), so two extensions cannot collide and no built-in is shadowed.
 */
export const ExtensionNamePattern = '^[a-z][a-z0-9_-]*:[a-z][a-z0-9_-]*$';

/** A built-in name or an extension's (`layout`, `heights`, `content`): the built-ins' enum or a namespaced name. */
const builtInOr = (names: readonly string[], doc: string) =>
  ({
    type: 'union',
    of: [
      { type: 'string', enum: names },
      { type: 'string', pattern: ExtensionNamePattern },
    ],
    doc,
  }) as const satisfies ParamType & { doc: string };

/**
 * A channel slot (_spec/README.md "Themes", Forms): a string (the `uniform` form: a color, a token, a word), a number
 * (a uniform size or opacity), or a node of `seq.theme`, read and printed in place.
 */
const channel = (doc: string) =>
  ({
    type: 'union',
    of: [{ type: 'string' }, { type: 'number' }, { type: 'node', language: ThemeLanguageId }],
    doc,
  }) as const satisfies ParamType & { doc: string };

export const SeqLanguage = defineLanguage({
  id: 'seq',
  version: SeqLanguageVersion,
  root: 'spec',
  doc: 'SeqViewSpec: what data there is (sequences, annotations, the maps between them) and how to show it (sections of tracks of layers).',
  kinds: {
    spec: {
      doc: 'The root of a spec.',
      parents: [],
      params: {
        version: {
          type: 'string',
          pattern: '^\\d+\\.\\d+\\.\\d+$',
          required: true,
          doc: 'The language version the spec is written in, `major.minor.patch`.',
        },
        title: { type: 'string', doc: 'A title for people.' },
        extensions: {
          type: 'map',
          of: { type: 'string', pattern: '^\\d+\\.\\d+(\\.\\d+)?$' },
          doc: "The extensions it names, by id, each the least version written for (`{ rondo: '1.2' }`).",
        },
      },
    },
    data: {
      doc: 'The container of the data of a spec; its id is `data`.',
      parents: ['spec'],
      ids: 'required',
    },
    views: {
      doc: 'The container of the views of a spec; its id is `views`.',
      parents: ['spec'],
      ids: 'required',
    },
    refs: {
      doc: 'The container of the user-defined refs of a spec (bookmarks, regions, sets); its id is `refs`.',
      parents: ['spec'],
      ids: 'required',
    },
    ref: {
      doc: 'A user-defined name for what an address names: `/refs/<id>` stands wherever an address does (a bookmark, a region, a set).',
      parents: ['refs'],
      ids: 'required',
      params: {
        to: {
          type: 'union',
          of: [
            { type: 'ref', doc: 'The address it names: a node, a node data makes, positions, items or a rect.' },
            {
              type: 'list',
              of: { type: 'ref', doc: 'An address of the set.' },
              'min-length': 1,
              doc: 'A set: what each address names, together.',
            },
          ],
          required: true,
          doc: 'What the ref names: an address, or a set of addresses (one type of result).',
        },
        label: { type: 'string', doc: 'A label for people (its id when absent).' },
        color: { type: 'string', doc: 'A color it is shown with: a hex literal or `token:<name>`.' },
      },
    },
    snapshots: {
      doc: 'The container of the story steps of a spec; its id is `snapshots`.',
      parents: ['spec'],
      ids: 'required',
    },
    snapshot: {
      doc: 'A story step (`/snapshots/<id>`): a full state, or a stored patch with no base shown over the spec (or over the step it starts `from`); entering it is navigation, never an edit.',
      parents: ['snapshots'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'A label for people (its id when absent).' },
        patch: {
          type: 'node',
          language: 'seq.patch',
          doc: 'What the step changes: a patch of this language with no `base`, its addresses refs of the spec (one of `patch` and `state`).',
        },
        state: {
          type: 'node',
          language: 'seq',
          kinds: ['spec'],
          doc: 'The step as a full state: a spec of its own, its refs resolved inside it; it does not follow edits of the spec (one of `patch` and `state`).',
        },
        from: {
          type: 'ref',
          kinds: ['snapshot'],
          doc: 'The step a patch applies over (a delta), instead of the spec.',
        },
        keep: {
          type: 'string',
          enum: ['full', 'spec'],
          doc: 'How the story’s encoder stores the step whatever is smaller: `full` a frozen full state, `spec` a patch that follows the spec.',
        },
        view: {
          type: 'record',
          fields: {
            view: { type: 'ref', kinds: ['view'], doc: 'The view shown (the first when absent).' },
            sections: {
              type: 'map',
              of: {
                type: 'record',
                fields: {
                  start: { type: 'number', required: true, doc: 'The first position shown (a text’s anchor).' },
                  end: { type: 'number', required: true, doc: 'The end of the window, a position.' },
                  top: { type: 'number', min: 0, doc: 'The vertical scroll, px.' },
                  scale: { type: 'number', min: 0, doc: 'px a position (a text’s letter size).' },
                  bands: {
                    type: 'map',
                    of: {
                      type: 'record',
                      fields: {
                        bandTop: { type: 'number', min: 0, doc: 'The band’s top row, in display rows.' },
                        rowScale: { type: 'number', min: 0, doc: 'px a row.' },
                      },
                    },
                    doc: 'A band that scrolls its own rows, by track address.',
                  },
                  y: {
                    type: 'record',
                    fields: {
                      start: { type: 'number', required: true, doc: 'Its first position.' },
                      end: { type: 'number', required: true, doc: 'Its end, a position.' },
                    },
                    doc: 'A grid’s second axis, as positions.',
                  },
                },
              },
              doc: 'Each section’s viewport by section address, as positions (never offsets).',
            },
          },
          doc: 'What the step shows: the view and each section’s viewport; a story names it on purpose, as data of the step.',
        },
        duration: { type: 'number', min: 0, doc: 'The transition into the step, ms (`story.transition.duration`).' },
        easing: {
          type: 'string',
          enum: ['linear', 'ease-in', 'ease-out', 'ease-in-out'],
          doc: 'How the transition’s time runs (`story.transition.easing`).',
        },
        linger: { type: 'number', min: 0, doc: 'How long the step is shown, ms (`story.linger`).' },
      },
    },
    source: {
      doc: 'Where data is read from: a file or service by uri, or a value inline.',
      parents: ['data'],
      ids: 'required',
      params: {
        format: {
          type: 'string',
          pattern: '^[a-z][a-z0-9_-]*$',
          required: true,
          doc: 'The format the source is read as (`fasta`, `bed`, `bigwig`).',
        },
        uri: { type: 'string', doc: 'Where it is: a url or a path.' },
        ref: {
          type: 'string',
          doc: 'A host’s name for the data (an edited or predicted model with no file), read by a host’s adapter only: beside a `uri`, which a reader without the host opens; alone, the source is the host’s (not portable, said once).',
        },
        index: { type: 'string', doc: 'Where its index is, for formats read by region.' },
        value: { type: 'any', doc: 'The data inline, for small values (instead of `uri`).' },
        fields: {
          type: 'map',
          of: { type: 'string', doc: 'The column name or the dotted path the source gives the field.' },
          doc: 'Which of the columns of the source is which field: a role (`target`, `start`, `end`, `position`, `strand`, `id`, `value`) or any other field of the ranges, by its column or its path in a record (`location.start.value`). Absent: each role from the member of its name, other scalar members fields under their names.',
        },
        assembly: { type: 'string', doc: 'The assembly a genomic source is on (`GRCh38`).' },
        records: {
          type: 'string',
          doc: 'The dotted path from the document’s root to its records (JSON, MessagePack): an array of records or an object of columns. Absent: the root.',
        },
        shape: {
          type: 'string',
          enum: ['records', 'matrix'],
          doc: 'What the source serves: `records` (absent), or a `matrix`, an array of arrays (JSON, MessagePack) found by `records` that a grid reads whole as its values, `null` missing.',
        },
        letters: {
          type: 'string',
          doc: 'The dotted path from the document’s root to the letters (a string) of the source’s own system. Absent: the source serves no letters.',
        },
        system: {
          type: 'string',
          doc: 'The name of the source’s own system: the rows with no `target`. Absent: the source node’s id.',
        },
        base: {
          type: 'integer',
          min: 0,
          max: 1,
          doc: 'The number of the file’s first position, `1` or `0`. Absent: `1` for a file, `0` for an inline `value`.',
        },
        ends: {
          type: 'string',
          enum: ['inclusive', 'exclusive'],
          doc: 'Whether a row’s `end` is its last position or the one after. Absent: `inclusive` for a file, `exclusive` for an inline `value`.',
        },
        types: {
          type: 'map',
          of: { type: 'string', enum: ['int32', 'float32', 'float64', 'category', 'string'] },
          doc: 'A type a field that is not a role: absent, a text cell is a `category`, a JSON number `float64`, a boolean a `category`.',
        },
        attribution: { type: 'string', doc: 'The data’s attribution line.' },
        licence: { type: 'string', doc: 'The data’s licence id (`CC-BY-4.0`).' },
      },
    },
    sequence: {
      doc: 'A sequence: letters of an alphabet, inline or from a source; a coordinate system.',
      parents: ['data'],
      ids: 'required',
      params: {
        alphabet: {
          type: 'string',
          enum: ['protein', 'dna', 'rna'],
          required: true,
          doc: 'The alphabet of its letters.',
        },
        letters: { type: 'string', doc: 'The letters inline (instead of `source`).' },
        source: { type: 'ref', kinds: ['source'], doc: 'The source that holds the letters (instead of `letters`).' },
        cyclic: { type: 'boolean', default: false, doc: 'Whether the sequence wraps around.' },
        name: {
          type: 'string',
          doc: 'The name its source knows it by (a FASTA record, a chromosome); the node’s id when absent.',
        },
        modified: ModifiedParam,
      },
    },
    numbering: {
      doc: 'A numbering of its sequence’s positions (its id the numbering’s name, `:numbering kabat`): a label a position, inline or from a source’s column. Labels need no order (IMGT), carry insertion codes (`100A`, `111.1`) and may be absent at a position.',
      parents: ['sequence', 'unit'],
      ids: 'required',
      params: {
        method: {
          type: 'string',
          enum: ['labels', 'offset', 'segments'],
          required: true,
          doc: 'How the labels are given: `labels`, one a position (inline or from a source); `offset`, a position plus `by`; `segments`, runs of consecutive numbers.',
        },
        by: {
          type: 'integer',
          doc: 'For `offset`: what is added to a position to give its label (1: numbered from 1).',
        },
        segments: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each run.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'The position after the last.' },
            first: { type: 'column', of: 'int32', required: true, doc: 'The number of its first position.' },
            insertion: {
              type: 'column',
              of: 'category',
              missing: true,
              doc: 'An insertion code after each run’s numbers (`52A`), for a run of one position.',
            },
          },
          doc: 'For `segments`: runs of positions numbered `first`, `first + 1`, …, not overlapping; a position in no run has no label.',
        },
        labels: {
          type: 'column',
          of: 'category',
          missing: true,
          doc: 'The labels inline, one a position of the sequence, in index order; `null` is no label there (instead of `source`).',
        },
        source: {
          type: 'ref',
          kinds: ['source'],
          doc: 'The source whose rows hold the labels, a row a position (its `position`, or its `start`), on the sequence (instead of `labels`).',
        },
        field: {
          type: 'string',
          doc: 'With `source`: the field of its rows that holds the label (`label` when absent).',
        },
        scheme: {
          type: 'string',
          pattern: '^[a-z][a-z0-9_-]*$',
          doc: 'The numbering scheme the labels follow, an open lowercase word (`^[a-z][a-z0-9_-]*$`: `kabat`, `imgt`, `chothia`, `uniprot`, `auth`, a lab’s own): what a reader says it is; it does not change the labels.',
        },
      },
    },
    ranges: {
      doc: 'Ranges on a sequence, 0-based and half-open: an annotation. On a cyclic sequence a range may start after it ends: it wraps the origin.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: {
          type: 'ref',
          kinds: ColumnsTargetKinds,
          required: true,
          doc: 'The system the ranges are on: a sequence, a molecule’s residues (every unit over it too), a unit, an alignment’s columns.',
        },
        rows: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each range.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'The position after the last, of each range.' },
            id: { type: 'column', of: 'category', doc: 'An identifier of each range.' },
            label: { type: 'list', of: { type: 'string', doc: 'A label.' }, doc: 'A label of each range.' },
            category: { type: 'column', of: 'category', doc: 'A category of each range (a theme colors by it).' },
            strand: { type: 'column', of: 'category', doc: 'The strand of each range: `+` or `-`.' },
          },
          doc: 'The ranges inline, a column for each field (instead of `source`).',
        },
        source: { type: 'ref', kinds: ['source'], doc: 'The source that holds the ranges (instead of `rows`).' },
        map: {
          type: 'string',
          enum: ['codon'],
          doc: '`codon`: each item (its rows, grouped by `id`) is a coding sequence on the target, in genomic order with its `strand`: the item gets a coordinate system of its residues and a codon map onto the target, made when something names it.',
        },
      },
    },
    values: {
      doc: 'A number per position of a sequence: a signal.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: {
          type: 'ref',
          kinds: SystemKinds,
          required: true,
          doc: 'The system the values are on: a sequence, a molecule’s residues, a unit.',
        },
        values: {
          type: 'column',
          of: 'float32',
          doc: 'One value per position, inline (instead of `source`); NaN is missing.',
        },
        source: {
          type: 'ref',
          kinds: ['source'],
          doc: 'The source that holds the values, read by region (instead of `values`).',
        },
        field: {
          type: 'string',
          doc: 'With `source`: the column of its values read (`plddt`, `b_factor`); the source’s default value when absent.',
        },
      },
    },
    alignment: {
      doc: 'An alignment: rows (a coordinate system of its own, in data order) over columns (the section’s system), inline or from a source.',
      parents: ['data'],
      ids: 'required',
      params: {
        alphabet: {
          type: 'string',
          enum: ['protein', 'dna', 'rna'],
          required: true,
          doc: 'The alphabet of its letters.',
        },
        rows: {
          type: 'table',
          columns: {
            label: {
              type: 'list',
              of: { type: 'string', doc: 'A label.' },
              required: true,
              doc: 'Each row’s label (its key): `NAME/27-137` numbers its first residue 27 (EC-39).',
            },
            letters: {
              type: 'list',
              of: { type: 'string', doc: 'A row’s aligned letters.' },
              required: true,
              doc: 'Each row’s aligned letters, every row as wide; gaps `-` and `.`.',
            },
          },
          doc: 'The rows inline (instead of `source`).',
        },
        source: {
          type: 'ref',
          kinds: ['source'],
          doc: 'The source that holds the alignment (an aligned FASTA or Stockholm file, or one serving rows by region), instead of `rows`.',
        },
        reference: {
          type: 'string',
          doc: 'The reference row’s label: what `identity` and `matches-reference` compare with (none when absent).',
        },
      },
    },
    tree: {
      doc: 'A tree over an alignment’s rows (Newick): its leaves match rows by label, else by the label before its last `/` (EC-43).',
      parents: ['data'],
      ids: 'required',
      params: {
        of: { type: 'ref', kinds: ['alignment'], required: true, doc: 'The alignment whose rows it orders.' },
        newick: { type: 'string', required: true, doc: 'The tree as Newick text.' },
      },
    },
    pairs: {
      doc: 'Pairs of locations that belong together (base pairs, contacts, links): each row two ends, a range on `a` and one on `b` (the same system when `b` is absent); held whole. As a child of a complex or a molecule, its links: it names no `a` or `b` (its parent is both), a complex’s rows each end’s unit.',
      parents: ['data', 'complex', 'molecule'],
      ids: 'required',
      params: {
        a: {
          type: 'ref',
          kinds: ColumnsTargetKinds,
          doc: 'The system each pair’s first end is on (an alignment: its columns, a family’s `SS_cons`); required under `data`, none under a complex or a molecule.',
        },
        b: {
          type: 'ref',
          kinds: ColumnsTargetKinds,
          doc: 'The system each pair’s second end is on (`a` when absent: pairs within one system).',
        },
        rows: {
          type: 'table',
          columns: {
            a_start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each first end.' },
            a_end: { type: 'column', of: 'int32', required: true, doc: 'The position after its last.' },
            b_start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each second end.' },
            b_end: { type: 'column', of: 'int32', required: true, doc: 'The position after its last.' },
            id: { type: 'column', of: 'category', doc: 'An identifier of each pair.' },
            type: {
              type: 'column',
              of: 'category',
              doc: 'What each pair is, an open word (`canonical`, `wobble`, `non-canonical`, `contact`): its kind (`kind` is the tree’s own word).',
            },
            label: { type: 'list', of: { type: 'string', doc: 'A label.' }, doc: 'A label of each pair.' },
            a_strand: { type: 'column', of: 'category', doc: 'The strand of each first end: `+` or `-`.' },
            b_strand: { type: 'column', of: 'category', doc: 'The strand of each second end: `+` or `-`.' },
            ...LinkColumns,
          },
          required: true,
          doc: 'The pairs inline, a column for each field; an end never wraps. Any other category column rides along as a field of the pairs a data filter reads (`(item :tertiary yes)`).',
        },
      },
    },
    grid: {
      doc: 'A value per pair of positions (a base-pair probability matrix, a PAE): `x`-major, held whole.',
      parents: ['data'],
      ids: 'required',
      params: {
        x: { type: 'ref', kinds: SystemKinds, required: true, doc: 'The system of the grid’s first axis.' },
        y: {
          type: 'ref',
          kinds: SystemKinds,
          doc: 'The system of its second axis (`x` when absent).',
        },
        values: {
          type: 'column',
          of: 'float32',
          required: true,
          doc: 'One value per pair, `x`-major (the value at x i, y j is at i × |y| + j); NaN is missing.',
        },
      },
    },
    map: {
      doc: 'A map from the positions of one sequence to those of another, so data on one shows on the other.',
      parents: ['data'],
      ids: 'required',
      params: {
        from: { type: 'ref', kinds: SystemKinds, required: true, doc: 'The system the map starts on.' },
        to: {
          type: 'ref',
          kinds: ColumnsTargetKinds,
          required: true,
          doc: 'The system the map ends on (an alignment: its columns, a chain onto the row it is closest to).',
        },
        method: {
          type: 'string',
          enum: ['offset', 'codon', 'segments', 'structure'],
          required: true,
          doc: '`offset`: a position plus `by`; `codon`: a protein on its gene through `exons`; `segments`: runs of positions (`segments`); `structure`: a sequence onto a structure’s residues (a unit or the entity’s molecule), SIFTS in abstract form: `segments` with each one’s `provenance`, and the residues’ `status`.',
        },
        by: { type: 'integer', doc: 'For `offset`: what is added to a position of `from`.' },
        segments: {
          type: 'table',
          columns: {
            from: { type: 'column', of: 'int32', required: true, doc: 'The first position of each run on `from`.' },
            to: { type: 'column', of: 'int32', required: true, doc: 'Its image on `to`.' },
            length: { type: 'column', of: 'int32', required: true, doc: 'The positions the run holds (at least 1).' },
            strand: {
              type: 'column',
              of: 'category',
              doc: 'Which way the run reads on `to`: `+` (when absent: `from + i` is `to + i`) or `-` (`to - i`).',
            },
            provenance: {
              type: 'column',
              of: 'category',
              missing: true,
              doc: 'For `structure`: where each run comes from, an open word (`sifts`, `alignment`, `author`).',
            },
          },
          doc: 'For `segments` and `structure`: the runs, not overlapping on `from`; what no run holds maps nowhere.',
        },
        status: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'The first residue on `to`.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'The residue after the last.' },
            status: {
              type: 'column',
              of: 'category',
              required: true,
              doc: '`conflict` (paired, but not the same residue: an engineered mutation, a microheterogeneity) or `unobserved` (modeled by no atoms; a unit’s `observed` wins).',
            },
          },
          doc: 'For `structure`: the residues of `to` that are not plainly paired; the rest of a run is paired, what no run holds unpaired (an expression tag, a linker).',
        },
        priority: {
          type: 'integer',
          doc: 'Its rank among maps joining the same systems (0 when absent): of the shortest paths that map differently, the one whose maps’ priorities, sorted weakest first and compared in turn, rank highest is taken; a full tie stays ambiguous.',
        },
        exons: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'The first position of each exon on `to`.' },
            end: { type: 'column', of: 'int32', required: true, doc: 'The position after the last, of each exon.' },
          },
          doc: 'For `codon`: the exons on `to`, in genomic order and not overlapping; `strand` says which way the protein runs. On a cyclic `to`, the last exon may end before it starts: it crosses the origin.',
        },
        strand: { type: 'string', enum: ['+', '-'], doc: 'For `codon`: the strand of the gene on `to`.' },
        table: {
          type: 'integer',
          min: 1,
          doc: 'For `codon`: the NCBI genetic code the codons are read with; 1 (the standard code) when absent.',
        },
      },
    },
    ...MoleculeKinds,
    view: {
      doc: 'A view: sections of tracks of layers, laid out as `layout` says.',
      parents: ['views'],
      ids: 'required',
      params: {
        layout: {
          ...builtInOr(
            ['tracks', 'text', 'rings', 'grid'],
            'How the sections are laid out (a section’s own `layout` before it): `tracks` (rows along one axis), `text` (lines of letters), `rings` (concentric rings, the first track outermost), `grid` (a plot of two axes).',
          ),
          default: 'tracks',
        },
        export: {
          type: 'record',
          fields: {
            legend: {
              type: 'string',
              enum: ['none', 'below', 'right'],
              doc: 'Legends in an export: `none` (when absent), in a band `below` the data or to its `right`.',
            },
            window: {
              type: 'boolean',
              doc: 'Whether an export draws a ring’s window (false when absent: an exported map is a figure).',
            },
            navigators: {
              type: 'string',
              enum: ['omit', 'include'],
              doc: 'Navigators in an export: `omit` (when absent: a window is session, not a figure), or `include`, stacked by `place` around the view, their windows as outlines.',
            },
          },
          doc: 'How the view exports (a figure says how it exports); an export call’s options override it.',
        },
      },
    },
    navigators: {
      doc: 'The container of a view’s navigators; its id is `navigators`.',
      parents: ['view'],
      ids: 'required',
    },
    navigator: {
      doc: 'A navigator: a part over a section’s axis, its rows or both, drawing the whole and the section’s window in it; a drag on the window moves the section. Placed by the host (`place` a hint).',
      parents: ['navigators'],
      ids: 'required',
      params: {
        of: {
          type: 'union',
          of: [
            { type: 'ref', kinds: ['section'] },
            { type: 'list', of: { type: 'ref', kinds: ['section'], doc: 'A section.' }, 'min-length': 1 },
          ],
          required: true,
          doc: 'The section, or sections, whose windows it draws and moves (a window each leading section; a follower’s window its leader’s carried along the follow).',
        },
        axes: {
          type: 'string',
          enum: ['positions', 'rows', 'both'],
          doc: 'The axes it covers: `positions` (when absent), an alignment’s `rows` (the window full width), or `both`.',
        },
        content: builtInOr(
          ['axis', 'tracks', 'bands', 'rows'],
          'What it draws of the whole: the `axis` as a hairline (when absent), the section’s `tracks` small, the `bands` of a chromosome (an ideogram), the section’s `rows` (summarized),.',
        ),
        tracks: {
          type: 'list',
          of: { type: 'ref', kinds: ['track'], doc: 'A track of a bound section.' },
          doc: 'For `tracks`: the tracks drawn small, in order; the section’s tracks but its rulers when absent.',
        },
        bands: {
          type: 'ref',
          kinds: ['ranges'],
          doc: 'For `bands`: the annotation of bands, a `stain` field each (UCSC’s cytoBandIdeo: `gneg`, `gpos25`…`gpos100`, `acen`, `gvar`, `stalk`).',
        },
        system: {
          type: 'ref',
          kinds: ['sequence'],
          doc: 'The system it draws when not the section’s: a chromosome under a section on a region of it, a gene under its protein. The section’s when absent.',
        },
        place: {
          type: 'string',
          enum: ['top', 'bottom', 'left', 'right'],
          doc: 'Where the host places it, a hint to the host and to an export: `top` when absent (`right` for `rows`).',
        },
        size: {
          type: 'union',
          of: [
            { type: 'number', min: 1 },
            { type: 'list', of: { type: 'number', min: 1 }, 'min-length': 2, 'max-length': 2 },
          ],
          doc: 'Its thickness in px (across its axis), or `[width, height]`: its form’s when absent (20 for an axis strip, 48 for tracks, 24 for bands, 96 for rows, 160 × 96 for both axes).',
        },
        hidden: { type: 'boolean', default: false, doc: 'Whether it is hidden: kept in the spec, not drawn.' },
      },
    },
    scale: {
      doc: 'A shared scale: layers whose themes name it share its domain and palette (an open end: the union over their data in view).',
      parents: ['views'],
      ids: 'required',
      params: ScaleParamsDescription,
    },
    section: {
      doc: 'A section of a view: a region of one coordinate system (or its `block` children), with its tracks.',
      parents: ['view'],
      ids: 'required',
      params: {
        of: {
          type: 'ref',
          kinds: ['sequence', 'alignment', 'molecule', 'unit', 'complex'],
          doc: 'The coordinate system the section shows (an alignment: its columns), or a complex (a block a unit, copies merged by `copies`, ordered by `classes`); required unless the section has `block` children.',
        },
        range: {
          type: 'list',
          of: { type: 'integer', min: 0, doc: 'A position.' },
          'min-length': 2,
          'max-length': 2,
          doc: 'The region shown, `[start, end)`; all of the system when absent.',
        },
        follow: {
          type: 'ref',
          kinds: ['section'],
          doc: 'A section whose region this one shows, mapped onto its own system.',
        },
        direction: {
          type: 'union',
          of: [
            {
              type: 'string',
              enum: ['genomic', 'reversed', 'follow'],
            },
            {
              type: 'record',
              fields: {
                along: {
                  type: 'ref',
                  kinds: ['sequence', 'map'],
                  required: true,
                  doc: 'The system (or a map from it) the section reads along.',
                },
              },
            },
          ],
          doc: 'Which way positions run left to right: `genomic` (up), `reversed` (down), `follow` (the way the leader reads, through the follow path), or `{ along: <system or map> }` (the way that system reads here). `follow` for a section that follows another, else `genomic`, when absent; never inferred from what the section draws.',
        },
        via: {
          type: 'union',
          of: [
            { type: 'ref', kinds: ['map'] },
            { type: 'list', of: { type: 'ref', kinds: ['map'], doc: 'A map.' }, 'min-length': 1 },
          ],
          doc: 'The map, or maps, this section’s paths take where other maps join the same two systems: its follow, its `direction: { along }` and what is projected onto it (hover, selection, the focus hull) leave out their rivals; other paths are unchanged.',
        },
        pieces: {
          type: 'union',
          of: [
            { type: 'string', enum: ['whole', 'coding'] },
            { type: 'ref', kinds: ['ranges', 'map'] },
            {
              type: 'list',
              of: { type: 'ref', kinds: ['ranges', 'map'], doc: 'A ranges node or a codon map.' },
              'min-length': 1,
            },
          ],
          doc: 'The positions shown, the rest hidden as breaks: `whole` (the system, or its `range`), `coding` (the coding runs of the codon map the section follows or reads along), a ranges node or a codon map, or a list of them (their union); intersected with `range`. `whole` when absent.',
        },
        collapsed: {
          type: 'union',
          of: [
            { type: 'string', enum: ['coding'] },
            { type: 'ref', kinds: ['ranges', 'map'] },
            {
              type: 'list',
              of: { type: 'ref', kinds: ['ranges', 'map'], doc: 'A ranges node or a codon map.' },
              'min-length': 1,
            },
          ],
          doc: 'The `pieces` the `section.collapse` command switches to from `whole` (`coding` when absent).',
        },
        pad: {
          type: 'integer',
          min: 0,
          doc: 'Positions kept each side of every piece (splice sites); overlapping pads merge. 0 when absent.',
        },
        layout: builtInOr(
          ['tracks', 'text', 'rings', 'grid'],
          'How the section is laid out: `tracks`, `text` (lines of letters, wrapped to the width), `rings` (its tracks as concentric rings round its whole axis), `grid` (a plot: `of` across, `of_y` down, its tracks’ layers over it); the view’s `layout` when absent.',
        ),
        of_y: {
          type: 'ref',
          kinds: ['sequence'],
          doc: 'Grid: the coordinate system down the plot (`of` when absent).',
        },
        follow_y: {
          type: 'union',
          of: [
            { type: 'string', enum: ['x', 'none'] },
            { type: 'ref', kinds: ['section'] },
          ],
          doc: 'Grid: what y shows: `x` (when absent: x’s window, so cells are square at aspect 1), a section whose region it shows (mapped onto `of_y`), or `none` (a window of its own).',
        },
        aspect: {
          type: 'union',
          of: [
            { type: 'number', min: 0 },
            { type: 'string', enum: ['none'] },
          ],
          doc: 'Grid: the plot’s width over its height (1 when absent: square cells), its height capped by the view’s room less a tenth; `none` leaves the height to the section’s tracks.',
        },
        group: {
          type: 'integer',
          min: 0,
          doc: 'Text: positions a group (a gap between groups; a line holds whole groups); 10 when absent, 0 for none.',
        },
        numbers: {
          type: 'string',
          enum: ['start', 'end', 'both', 'none'],
          doc: 'Text: each line’s first (`start`, left) or last (`end`, right) position’s number, or both, or none; `start` when absent.',
        },
        group_labels: {
          type: 'string',
          enum: ['above', 'inline', 'none'],
          doc: 'Text: a number a group, above its last position (`above`, when absent), before the group on its line (`inline`), or none.',
        },
        align: {
          type: 'string',
          enum: ['start', 'center'],
          doc: 'Text: lines at the start of the width (when absent) or centered in it.',
        },
        max_width: { type: 'number', min: 0, doc: 'Text: the most px a line takes; the whole width when absent.' },
        height: {
          type: 'number',
          min: 0,
          doc: 'Text: the section’s height in px (it scrolls inside); its content’s, up to the render config’s `text.maxHeight`, when absent. Rings: the rings’ diameter in px; the width’s, up to `rings.maxDiameter`, when absent.',
        },
        heights: {
          type: 'union',
          of: [
            { type: 'string', enum: ['content', 'uniform', 'fixed'] },
            { type: 'string', pattern: ExtensionNamePattern },
          ],
          doc: 'How the tracks with no `height` are sized: `content` (when absent: each as tall as its lanes over the whole data, in steps of a lane, tracks of one kind the tallest of them), `uniform` (every such track the tallest content height in the section), `fixed` (each its layer type’s fixed height), or an extension’s `<id>:<name>`. A track’s `height` overrides it.',
        },
        rotation: {
          type: 'number',
          doc: 'Rings: where the origin is, degrees clockwise from the top; 0 when absent.',
        },
        labels: {
          type: 'string',
          enum: ['outside', 'none'],
          doc: 'Rings: the items’ labels placed outside the rings with leaders (`outside`, when absent), or none.',
        },
        arrange: {
          type: 'string',
          enum: ['row', 'stack'],
          doc: 'Blocks in tracks: side by side along one axis (`row`, when absent) or a row of tracks each, aligned at their first position (`stack`).',
        },
        copies: {
          type: 'string',
          enum: ['merge', 'each'],
          doc: 'Blocks that are copies (the same system and range): one block labeled `×N` (`merge`, when absent) or a block each.',
        },
        classes: SectionClasses,
      },
    },
    block: {
      doc: 'A block of a section: a piece of it on a system of its own (a chain of a complex, a ligand); blocks follow one another in order.',
      parents: ['section'],
      ids: 'required',
      params: {
        of: {
          type: 'ref',
          kinds: SystemKinds,
          doc: 'The coordinate system the block shows: a sequence, a molecule, a unit of a complex. Absent on a block of a section of a complex whose id is a unit’s: it overrides that unit (its `form`, its `label`), the others still derived.',
        },
        range: {
          type: 'list',
          of: { type: 'integer', min: 0, doc: 'A position.' },
          'min-length': 2,
          'max-length': 2,
          doc: 'The region shown, `[start, end)` (over a cyclic origin when start is past end); all of the system when absent.',
        },
        label: { type: 'string', doc: 'The block’s label (its id when absent).' },
        form: BlockForm,
      },
    },
    track: {
      doc: 'A track of a section: a row of layers.',
      parents: ['section'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'The label of the track (its id when absent).' },
        height: { type: 'number', min: 0, doc: 'The height of the track in px.' },
        padding: {
          type: 'union',
          of: [
            { type: 'number', min: 0 },
            {
              type: 'record',
              fields: {
                top: { type: 'number', min: 0, doc: 'Px above the data.' },
                bottom: { type: 'number', min: 0, doc: 'Px below the data.' },
              },
            },
          ],
          doc: 'Px inside the track’s well around its data, vertical only: a number is top and bottom, `{ top, bottom }` one or both sides (the render config’s `trackPadding` when absent, 0). The well is the track’s height plus its padding; its layers lay out, hit and are emphasized in the band inside it (a ring’s padding is radial).',
        },
        hidden: {
          type: 'boolean',
          default: false,
          doc: 'Whether the track is hidden: kept in the spec, its row collapsed; its layers keep their own `hidden`.',
        },
        for_each: {
          type: 'record',
          fields: {
            source: {
              type: 'ref',
              kinds: ['source'],
              required: true,
              doc: 'The source whose dynamic nodes the template shows.',
            },
            node_kind: {
              type: 'string',
              enum: ['ranges', 'values'],
              required: true,
              doc: 'The kind of node a track is made for: an annotation, which its layers draw.',
            },
            depth: {
              type: 'integer',
              min: 0,
              doc: 'The nodes of the sources nested this many levels below `source` (0, its own, when absent).',
            },
            order: {
              type: 'string',
              pattern: '^-?[a-z_$][a-z0-9_]*$',
              doc: 'The order of the tracks: `key` (when absent) or a param of the nodes, descending with a leading `-`.',
            },
            filter: {
              type: 'map',
              of: { type: 'string', doc: 'The value the param must have, as written.' },
              doc: 'The params a node must have to get a track, each the value given.',
            },
          },
          doc: 'A template: a track for each dynamic node of a kind a source holds, in order, its layers drawing that node (a layer without `data` draws the node).',
        },
      },
    },
    layer: {
      doc: 'A layer of a track: data drawn as one type of mark.',
      parents: ['track'],
      ids: 'required',
      // an extension's type (`rondo:ribbons`) reads plain params beside these: kept, checked by the extension
      extensionParam: 'type',
      params: {
        type: {
          // a word, not an enum: a type this reader lacks (an extension's, a newer built-in) is drawn as its fallback and
          // said by the resolved view (`seq.view.layer.type`, `seq.view.layer.extension`), never refused
          type: 'string',
          pattern: '^([a-z][a-z0-9_-]*:)?[a-z][a-z0-9_-]*$',
          required: true,
          doc: 'What it draws (an extension’s `<id>:<name>` too): the letters of the system, ranges, values, the ruler of the section’s axis, secondary structure (ranges whose `kind` field says helix, strand or turn), variants (ranges of one position each), frames (ranges outlined across the tracks `across` names, titled by their labels), an alignment’s rows (`rows`) or its conservation logo (`logo`), pairs (`pairs`: arcs over a track, chords on a ring, dots on a grid), a grid of values (`grid`, in a grid section) or a chromosome’s stained bands (`bands`, an ideogram: ranges with a `stain` field, colored by the `giemsa` scheme, the centromere’s two `acen` bands as triangles); `secondary`, `variants` and `bands` draw `ranges` where a reader lacks them, `pairs` draws its ends as `ranges`.',
        },
        fallback: {
          type: 'string',
          pattern: '^[a-z][a-z0-9_-]*$',
          doc: 'The built-in type a reader without `type` draws instead, or `none`.',
        },
        data: {
          type: 'ref',
          kinds: ['sequence', 'ranges', 'values', 'alignment', 'pairs', 'grid'],
          doc: 'The data it draws (the letters of the section for `letters`, its axis for `ruler`; a `ranges` layer takes `pairs` too: each end a mark at its site).',
        },
        color: channel('The marks’ color: a color, a token, `light-dark(a, b)`, or a theme form.'),
        ink: channel('`letters`: the letters’ color (`auto`: follows its fill, 4.5:1), or a theme form.'),
        size: channel('The marks’ size in the layer’s unit (a share of the slot, or px), or a theme form.'),
        opacity: channel('The marks’ opacity, 0 to 1 (never a way to hide: a filter hides), or a theme form.'),
        line: channel(
          '`ranges`: how an item is drawn as a line: `solid`, `dashed`, `dotted`, `wavy`, `double`, `box` or `highlight`; blocks when absent. Any other value is drawn solid and said once (`seq.render.line-style`).',
        ),
        arrow: channel(
          '`ranges`: an arrowhead on the item’s end, start or both (`end`, `start`, `both`), on the piece holding that end, never on a cut one; `none` when absent.',
        ),
        shape: channel('`variants`: the head’s shape (`circle`, `square`, `triangle`, `diamond`).'),
        scale: {
          type: 'ref',
          kinds: ['scale'],
          doc: '`values`: the scale its y reads (shared with the layers naming it).',
        },
        domain: {
          type: 'list',
          of: { type: 'number', doc: 'An end.' },
          'min-length': 2,
          'max-length': 2,
          doc: '`values`: its y domain inline, `[a, b]` (open over the values in view when absent).',
        },
        hidden: { type: 'boolean', default: false, doc: 'Whether the layer is hidden: kept in the spec, not drawn.' },
        complement: {
          type: 'string',
          enum: ['auto', 'always', 'never'],
          doc: 'For `letters` of nucleotides: `auto` (when absent) shows a letter as its complement where it is read against its strand on screen (the axis reversed against its run); `always`, `never`. The stored letter is what search and status read.',
        },
        across: {
          type: 'list',
          of: { type: 'string', doc: 'A track’s id in the layer’s section.' },
          'min-length': 1,
          'max-length': 2,
          doc: 'For `frame`: the tracks its outlines cover, the first and the last (the bands between them too); its own track’s when absent.',
        },
        numbering: {
          type: 'string',
          doc: 'For `ruler`: the numbering it labels by (`ordinal`, `index`, or one the system declares); the section’s default numbering (`ordinal`) when absent.',
        },
        units: {
          type: 'string',
          enum: ['none', 'bases'],
          doc: 'For `ruler`: `none` (when absent: the number as its numbering labels it) or `bases` (bp, kb, Mb, chosen by the tick step: 7,668,402 at a base tick, 7.67 Mb at a megabase one); a numbering other than `ordinal` keeps its own labels.',
        },
        unnumbered: {
          type: 'string',
          enum: ['mark', 'fallback'],
          doc: 'For `ruler`: what stands where its numbering has no label (past its coverage: the constant domain under Kabat): `mark` (when absent) a muted dashed line, one muted “not numbered” note at the start of each run and no ticks; `fallback` the run labelled by `fallbackNumbering`, muted.',
        },
        fallbackNumbering: {
          type: 'string',
          doc: 'For `ruler` with `unnumbered: fallback`: the numbering that labels the unnumbered positions (`ordinal`, `index`, or one the system declares); `ordinal` when absent. One the system lacks is `seq.render.numbering.fallback` and `mark` is used.',
        },
        encoding: {
          type: 'string',
          enum: ['bars', 'line', 'area'],
          doc: 'For `values`: `bars` (when absent), `line` (through the cells’ centers, held over half cells at its ends and gaps, broken by missing values) or `area` (the line and the area under it).',
        },
        style: {
          type: 'string',
          enum: ['blocks'],
          doc: 'For `secondary`: the figure’s form; `blocks` (when absent): helices and strands as blocks, a strand ending in an arrowhead.',
        },
        tree: {
          type: 'ref',
          kinds: ['tree'],
          doc: 'For `rows`: the tree drawn beside the rows and read by `order: tree`.',
        },
        order: {
          type: 'union',
          of: [
            { type: 'string', enum: ['data', 'tree'] },
            {
              type: 'record',
              fields: {
                field: {
                  type: 'string',
                  required: true,
                  doc: 'A field of the rows table (`label`, `first`, an own field, a measure, a joined field).',
                },
                direction: { type: 'string', enum: ['ascending', 'descending'], doc: '`ascending` when absent.' },
              },
            },
          ],
          doc: 'For `rows`: the rows’ order: `data` (when absent), `tree` (the tree’s leaf order; clades fold only there), or a field of the rows table (missing last, ties in data order).',
        },
        collapsed: {
          type: 'list',
          of: {
            type: 'string',
            doc: 'A clade: its node’s label where unique in the tree and not a number, else `first..last` (its first and last leaves).',
          },
          doc: 'For `rows`: the clades folded when the view opens (session state after: `clades.keep` writes the folds shown here).',
        },
        pinned: {
          type: 'list',
          of: { type: 'string', doc: 'A row’s label.' },
          doc: 'For `rows`: rows drawn in a band at the top (and kept in their place in the order, marked).',
        },
        hidden_rows: {
          type: 'list',
          of: { type: 'string', doc: 'A row’s label.' },
          doc: 'For `rows`: rows not shown (a count chip at the header’s foot shows them again).',
        },
        height: {
          type: 'union',
          of: [
            { type: 'string', enum: ['fill', 'span'] },
            { type: 'integer', min: 1 },
          ],
          doc: 'For `rows`: the band’s height, `fill` (when absent: the section’s room left) or a number of rows. For `pairs`: an arc’s height, `span` (when absent: by its span against the widest in view, so nested arcs nest).',
        },
        side: {
          type: 'string',
          enum: ['above', 'below'],
          doc: 'For `pairs`: arcs rise from the band’s foot (`above`, when absent) or hang from its top (`below`); on a ring `below` reaches the hole (chords).',
        },
        reach: {
          type: 'number',
          min: 0,
          doc: 'For `pairs`: how far past the window an end out of view is taken, a share of the window (0.25 when absent), so an arc moves continuously on pan.',
        },
        crossing: {
          type: 'string',
          enum: ['same', 'split'],
          doc: 'For `pairs`: crossing pairs drawn crossing (`same`, when absent) or `split`: the largest nested set on `side`, the rest on the other side.',
        },
        across_lines: {
          type: 'boolean',
          doc: 'For `pairs` in text: arcs between lines too (false when absent: a line’s own pairs, the others’ ends marked).',
        },
        mirror: {
          type: 'boolean',
          doc: 'For `pairs` in a grid: each pair drawn at (a, b) and at (b, a) (false when absent: where its data puts it).',
        },
        triangle: {
          type: 'string',
          enum: ['both', 'upper', 'lower'],
          doc: 'For `grid` and `pairs` in a grid: the whole plot (`both`, when absent), or only above (`upper`: y before x) or below its diagonal.',
        },
        aggregate: {
          type: 'string',
          enum: ['mean', 'max', 'min'],
          doc: 'For `grid`: what a sample reads of its cells zoomed out (`mean` when absent; `min` keeps a PAE’s confident pairs).',
        },
        summary: {
          type: 'string',
          enum: ['consensus', 'bar'],
          doc: 'For `rows`: a folded clade’s cells: `consensus` (when absent: the majority letter, opaque by its share) or `bar` (the non-gap share).',
        },
        numbers: {
          type: 'string',
          enum: ['gutters', 'none'],
          doc: 'For `rows`: each row’s first and last residue numbers in gutters (`gutters`, when absent) or none.',
        },
        columns: {
          type: 'list',
          of: {
            type: 'record',
            fields: {
              field: {
                type: 'string',
                required: true,
                doc: 'A field of the rows table: the alignment’s own (Stockholm’s `#=GS` tags: `OS`, `DE`, `AC`; an inline table’s other columns), a measure (`length`, `gaps`, `identity` to the reference row), or a field of a records source whose system is the rows system, joined by the row’s key.',
              },
              side: {
                type: 'string',
                enum: ['left', 'right'],
                doc: '`left` (when absent: in the header, between the label and the first numbers) or `right` (after the last numbers, before a rows rail).',
              },
              as: {
                type: 'string',
                enum: ['text', 'number', 'bar', 'chip'],
                doc: '`text` (when absent), `number` (right-aligned), `bar` (the value’s share of the column’s range) or `chip` (a category’s color).',
              },
              width: {
                type: 'number',
                min: 1,
                doc: 'Its width, px (by `as` when absent: text 96, number 44, bar 56, chip 56).',
              },
              title: {
                type: 'string',
                doc: 'Its title above the rows (the field when absent): a click orders the rows by it.',
              },
            },
          },
          doc: 'For `rows`: data columns beside each row’s header, their titles above the rows; on a narrow view they yield before the label, right ones first, each down to a floor and then hidden (a count chip says how many).',
        },
        tree_width: {
          type: 'number',
          min: 0,
          doc: 'For `rows`: the tree’s width in px (120 when absent: the view’s label column widens by it; 0 hides it).',
        },
        logo_scale: {
          type: 'string',
          enum: ['bits', 'probability'],
          doc: 'For `logo`: letter heights by information (`bits`, when absent: log2(k) − H) or by `probability`.',
        },
        correction: {
          type: 'string',
          enum: ['none', 'small-sample'],
          doc: 'For `logo` in `bits`: the small-sample correction (`none` when absent: exact for given frequencies).',
        },
        gaps: {
          type: 'string',
          enum: ['ignore', 'count'],
          doc: 'For `logo`: gaps left out of the frequencies (`ignore`, when absent) or counted.',
        },
        code: {
          type: 'string',
          enum: ['auto', 'one', 'three'],
          doc: 'For `letters` of residues a codon wide: `auto` (when absent) three-letter codes where the narrowest part in view has room; `one`, `three`.',
        },
        modified: LettersModified,
      },
    },
  },
});

export type SeqLanguageType = typeof SeqLanguage;

/** The language of the spec's stored patches (a story step's `patch`): `seq.patch`, `repr/patch`'s over this one. */
export const SeqPatchLanguage: Language = patchLanguage(SeqLanguage);

/**
 * The languages the spec's slots name, by id: its stored patches' (_spec/README.md, Language, Story steps) and the
 * theme sublanguage (_spec/README.md "Themes").
 */
export function seqLanguages(id: string): Language | undefined {
  // a full-state story step is a spec of its own (`snapshot.state`)
  if (id === SeqLanguage.id) return SeqLanguage;
  return id === SeqPatchLanguage.id ? SeqPatchLanguage : id === ThemeLanguageId ? ThemeLanguage : undefined;
}

/** What validation, normalization and `applyPatch` take to read a spec's slots (`validate` in `ApplyOptions`). */
export const SeqValidateOptions: ValidateOptions = { languages: seqLanguages };
export type SeqKinds = SeqLanguageType['kinds'];
