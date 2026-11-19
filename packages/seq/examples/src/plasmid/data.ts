/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { printAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { deriveFeatures, findSites, PlasmidFiles, PlasmidName, rotatedAnnotations } from './derive.js';
import { PublishedLetters } from './record.js';

export { PlasmidFiles, PlasmidLength, PlasmidName, RotatedOrigin } from './derive.js';

// The Plasmid dataset (examples.md, Plasmid): pUC19c, INSDC L09137.2 (recorded as a `seq` file, its letters also a
// module: `record.ts`), its features found in the sequence and its restriction sites computed on it, in code
// (`derive.ts`), the rows written inline; the views read the letters from the `seq` files. Two instances: the published
// one, and the rotated one whose origin is inside bla, so bla and a generated primer site cross it (its letters a `seq`
// file of their own). The ring's scale variant is generated here, from a seed. The rotated variant's features are also
// served as GFF3 (`fixtures/plasmid/`, the Files dataset's view reads it).

export type PlasmidVariant = 'published' | 'rotated';

type Builder = ReturnType<typeof createSeqSpec>;

/** A feature or a site as a ranges row: a feature across the origin ends before it starts (it wraps). */
const endOf = (f: { start: number; length: number }, length: number) => {
  const end = f.start + f.length;
  return end > length ? end - length : end;
};

/**
 * The sequence and the ranges of an instance: each strand's features, the sites, and the sites' address (all of them,
 * or the enzymes named). The letters are read from the instance's `seq` file; the annotations are derived from the
 * record's letters here, at once: the features are found in the sequence, the sites computed on it, never read.
 */
export function plasmidData(b: Builder, variant: PlasmidVariant = 'rotated') {
  const rotated = variant === 'rotated';
  const letters = PublishedLetters;
  const { features, sites } = rotated
    ? rotatedAnnotations(letters)
    : { features: deriveFeatures(letters), sites: findSites(letters) };
  const length = letters.length;
  const file = b.data.source('record', { format: 'seq', uri: rotated ? PlasmidFiles.rotated : PlasmidFiles.published });
  const seq = b.data.sequence('plasmid', { alphabet: 'dna', source: file, cyclic: true, name: PlasmidName });
  const strand = (s: '+' | '-') => {
    const f = features.filter((x) => x.strand === s);
    return b.data.ranges(s === '+' ? 'plus' : 'minus', {
      target: seq,
      rows: {
        start: f.map((x) => x.start),
        end: f.map((x) => endOf(x, length)),
        id: f.map((x) => x.id),
        label: f.map((x) => x.label),
        category: f.map((x) => x.type),
        strand: f.map(() => s),
      },
    });
  };
  const siteRanges = b.data.ranges('sites', {
    target: seq,
    rows: {
      // a site over its recognition sequence (one across the origin ends before it starts; two enzymes' can overlap:
      // lanes)
      start: sites.map((x) => x.start),
      end: sites.map((x) => endOf(x, length)),
      id: sites.map((x) => x.enzyme),
      label: sites.map((x) => x.enzyme),
    },
  });
  return {
    seq,
    plus: strand('+'),
    minus: strand('-'),
    /** The sites' address, those of the enzymes named only (a filter on the layer's data) or all. */
    sites: (enzymes?: readonly string[]) =>
      enzymes ? `${siteRanges}${printAddress([{ kind: 'item', params: { id: [...enzymes] } }])}` : siteRanges,
    /** The enzymes with a site on the instance. */
    enzymes: [...new Set(sites.map((x) => x.enzyme))],
  };
}

/** The ring's scale variant: 200 kb cyclic, 2,000 features on both strands, a tenth across the origin. */
export const ScaleLength = 200_000;
export const ScaleFeatures = 2000;

export interface ScaleData {
  readonly letters: string;
  readonly features: {
    readonly start: number[];
    readonly end: number[];
    readonly id: string[];
    readonly label: string[];
    readonly category: string[];
    readonly strand: ('+' | '-')[];
  };
}

/** Generated from a seed on core's `PCG`, never stored: the letters and the features (a ring's rows, `end` past `start`). */
export function scaleData(seed = 200): ScaleData {
  const r = new PCG(seed);
  let letters = '';
  for (let i = 0; i < ScaleLength; i++) letters += 'ACGT'[r.int() & 3];
  const f: { -readonly [K in keyof ScaleData['features']]: ScaleData['features'][K] } = {
    start: [],
    end: [],
    id: [],
    label: [],
    category: [],
    strand: [],
  };
  const types = ['CDS', 'misc_feature', 'rep_origin', 'primer_bind'];
  const crossing = ScaleFeatures / 10;
  for (let i = 0; i < ScaleFeatures; i++) {
    const length = 100 + Math.floor(r.float() * 1400);
    // a tenth of them start before the origin and end after it (a range that starts after it ends wraps)
    const start =
      i < crossing
        ? ScaleLength - 1 - Math.floor(r.float() * (length - 1))
        : Math.floor(r.float() * (ScaleLength - length));
    const end = i < crossing ? (start + length) % ScaleLength : start + length;
    f.start.push(start);
    f.end.push(end);
    f.id.push(`f${i}`);
    f.label.push(`feature ${i}`);
    f.category.push(types[i % types.length]);
    f.strand.push(r.int() & 1 ? '+' : '-');
  }
  return { letters, features: f };
}

/** The scale variant's ranges of one strand, inline. */
export function scaleStrand(b: Builder, seq: string, data: ScaleData, strand: '+' | '-', id: string) {
  const k = data.features.strand.map((s, i) => (s === strand ? i : -1)).filter((i) => i >= 0);
  const pick = <T>(a: readonly T[]) => k.map((i) => a[i]);
  return b.data.ranges(id, {
    target: seq,
    rows: {
      start: pick(data.features.start),
      end: pick(data.features.end),
      id: pick(data.features.id),
      label: pick(data.features.label),
      category: pick(data.features.category),
      strand: pick(data.features.strand),
    },
  });
}
