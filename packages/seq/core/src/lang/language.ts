/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqViewSpec, described as data on `@molstar/core/repr` (_spec/README.md, Language): the first slice's kinds. What
// the description cannot say (a sequence has letters or a source, a codon map links a protein to a gene) is the checks
// of `check.ts`, in the `seq.lang.*` codes.

import { patchLanguage } from '@molstar/core/repr/patch';
import { defineLanguage, type Language, type ValidateOptions } from '@molstar/core/repr/schema';

/**
 * The version of the language this description is (the package's `LanguageVersion` is the same). Pre-release (`0.x`):
 * a change moves it with no migration (_spec/README.md, Versions).
 */
export const SeqLanguageVersion = '0.4.0';

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
      doc: 'A story step: a stored patch with no base, shown over whatever the spec is (`/snapshots/<id>`); entering it is navigation, never an edit.',
      parents: ['snapshots'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'A label for people (its id when absent).' },
        patch: {
          type: 'node',
          language: 'seq.patch',
          required: true,
          doc: 'What the step changes: a patch of this language with no `base`, its addresses refs of the spec.',
        },
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
        index: { type: 'string', doc: 'Where its index is, for formats read by region.' },
        value: { type: 'any', doc: 'The data inline, for small values (instead of `uri`).' },
        fields: {
          type: 'map',
          of: { type: 'string', doc: 'The name the source gives the field.' },
          doc: 'Which of the columns of the source is which field.',
        },
        assembly: { type: 'string', doc: 'The assembly a genomic source is on (`GRCh38`).' },
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
      },
    },
    ranges: {
      doc: 'Ranges on a sequence, 0-based and half-open: an annotation. On a cyclic sequence a range may start after it ends: it wraps the origin.',
      parents: ['data'],
      ids: 'required',
      params: {
        target: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the ranges are on.' },
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
        target: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the values are on.' },
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
      },
    },
    map: {
      doc: 'A map from the positions of one sequence to those of another, so data on one shows on the other.',
      parents: ['data'],
      ids: 'required',
      params: {
        from: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the map starts on.' },
        to: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The sequence the map ends on.' },
        method: {
          type: 'string',
          enum: ['offset', 'codon'],
          required: true,
          doc: '`offset`: a position plus `by`; `codon`: a protein on its gene through `exons`.',
        },
        by: { type: 'integer', doc: 'For `offset`: what is added to a position of `from`.' },
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
    view: {
      doc: 'A view: sections of tracks of layers, laid out as `layout` says.',
      parents: ['views'],
      ids: 'required',
      params: {
        layout: { type: 'string', enum: ['tracks'], default: 'tracks', doc: 'How the sections are laid out.' },
      },
    },
    section: {
      doc: 'A section of a view: a region of one coordinate system, with its tracks.',
      parents: ['view'],
      ids: 'required',
      params: {
        of: { type: 'ref', kinds: ['sequence'], required: true, doc: 'The coordinate system the section shows.' },
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
      },
    },
    track: {
      doc: 'A track of a section: a row of layers.',
      parents: ['section'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'The label of the track (its id when absent).' },
        height: { type: 'number', min: 0, doc: 'The height of the track in px.' },
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
      params: {
        type: {
          type: 'string',
          enum: ['letters', 'ranges', 'values', 'ruler', 'secondary', 'variants'],
          required: true,
          doc: 'What it draws: the letters of the system, ranges, values, the ruler of the section’s axis, secondary structure (ranges whose `kind` field says helix, strand or turn) or variants (ranges of one position each); `secondary` and `variants` draw `ranges` where a reader lacks them.',
        },
        data: {
          type: 'ref',
          kinds: ['sequence', 'ranges', 'values'],
          doc: 'The data it draws (the letters of the section for `letters`, its axis for `ruler`).',
        },
        color: { type: 'any', doc: 'A theme form that colors the marks.' },
        hidden: { type: 'boolean', default: false, doc: 'Whether the layer is hidden: kept in the spec, not drawn.' },
        complement: {
          type: 'string',
          enum: ['auto', 'always', 'never'],
          doc: 'For `letters` of nucleotides: `auto` (when absent) shows a letter as its complement where it is read against its strand on screen (the axis reversed against its run); `always`, `never`. The stored letter is what search and status read.',
        },
        numbering: {
          type: 'string',
          doc: 'For `ruler`: the numbering it labels by (`ordinal`, `index`, or one the system declares); the section’s default numbering (`ordinal`) when absent.',
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
        line: {
          type: 'string',
          doc: 'For `ranges`: how an item is drawn as a line: `solid`, `dashed`, `dotted`, `wavy`, `double`, `box` or `highlight`; blocks when absent. Any other value is drawn solid and said once (`seq.render.line-style`).',
        },
        arrow: {
          type: 'string',
          enum: ['none', 'end', 'start', 'both'],
          doc: 'For `ranges`: an arrowhead on the item’s end, start or both, on the piece holding that end, never on a cut one; `none` when absent.',
        },
        code: {
          type: 'string',
          enum: ['auto', 'one', 'three'],
          doc: 'For `letters` of residues a codon wide: `auto` (when absent) three-letter codes where the narrowest part in view has room; `one`, `three`.',
        },
      },
    },
  },
});

export type SeqLanguageType = typeof SeqLanguage;

/** The language of the spec's stored patches (a story step's `patch`): `seq.patch`, `repr/patch`'s over this one. */
export const SeqPatchLanguage: Language = patchLanguage(SeqLanguage);

/** The languages the spec's slots name, by id: its stored patches' (_spec/README.md, Language, Story steps). */
export function seqLanguages(id: string): Language | undefined {
  return id === SeqPatchLanguage.id ? SeqPatchLanguage : undefined;
}

/** What validation, normalization and `applyPatch` take to read a spec's slots (`validate` in `ApplyOptions`). */
export const SeqValidateOptions: ValidateOptions = { languages: seqLanguages };
export type SeqKinds = SeqLanguageType['kinds'];
