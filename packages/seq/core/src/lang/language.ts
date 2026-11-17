/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqViewSpec, described as data on `@molstar/core/repr` (_spec/README.md, Language): the first slice's kinds. What
// the description cannot say (a sequence has letters or a source, a codon map links a protein to a gene) is the checks
// of `check.ts`, in the `seq.lang.*` codes.

import { defineLanguage } from '@molstar/core/repr/schema';

/** The version of the language this description is (the package's `LanguageVersion` is the same). */
export const SeqLanguageVersion = '0.0.0';

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
      },
    },
    ranges: {
      doc: 'Ranges on a sequence, 0-based and half-open: an annotation.',
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
            strand: { type: 'column', of: 'category', doc: 'The strand of each range, `+` or `-`.' },
          },
          doc: 'The ranges inline, a column for each field (instead of `source`).',
        },
        source: { type: 'ref', kinds: ['source'], doc: 'The source that holds the ranges (instead of `rows`).' },
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
          doc: 'For `codon`: the exons on `to`, in the order of the protein.',
        },
        strand: { type: 'string', enum: ['+', '-'], doc: 'For `codon`: the strand of the gene on `to`.' },
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
      },
    },
    track: {
      doc: 'A track of a section: a row of layers.',
      parents: ['section'],
      ids: 'required',
      params: {
        label: { type: 'string', doc: 'The label of the track (its id when absent).' },
        height: { type: 'number', min: 0, unit: 'px', doc: 'The height of the track in pixels.' },
      },
    },
    layer: {
      doc: 'A layer of a track: data drawn as one type of mark.',
      parents: ['track'],
      params: {
        type: {
          type: 'string',
          enum: ['letters', 'ranges', 'values'],
          required: true,
          doc: 'What it draws: the letters of the system, ranges, or values.',
        },
        data: {
          type: 'ref',
          kinds: ['sequence', 'ranges', 'values'],
          doc: 'The data it draws (the letters of the section for `letters`).',
        },
        color: { type: 'any', doc: 'A theme form that colors the marks.' },
      },
    },
  },
});

export type SeqLanguageType = typeof SeqLanguage;
export type SeqKinds = SeqLanguageType['kinds'];
