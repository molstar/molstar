/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, type Node, type Value } from '@molstar/core/repr/tree';
import { type SeqData, seqToValue } from '@molstar/seqx-core/data';
import {
  createSeqSpec,
  ideogramNavigator,
  miniatureNavigator,
  type NavigatorParams,
  stripNavigator,
} from '@molstar/seqx-core/lang';
import * as globins from './globins/family.js';
import * as geneByRegion from './p53/gene-by-region.js';
import * as proteinOnGene from './p53/protein-on-gene.js';
import { RegionSequence } from './p53/region.js';
import * as tracks from './p53/tracks.js';
import * as plasmidRing from './plasmid/ring.js';

// The examples' navigators: each a view of the examples with navigators added (a node under its view's
// `navigators`), and for the ideogram the chromosome and its bands besides.

/** chr17's length in GRCh38, positions. */
export const Chr17Length = 83_257_441;

/**
 * chr17's ideogram bands, GRCh38: UCSC's cytoBandIdeo (`hg38`, fetched from its public API, 2026-11-19), as given,
 * half-open, 0-based.
 */
export const Chr17Bands = {
  start: [
    0, 3400000, 6500000, 10800000, 16100000, 22700000, 25100000, 27400000, 33500000, 39800000, 40200000, 42800000,
    46800000, 49300000, 52100000, 59500000, 60200000, 63100000, 64600000, 66200000, 69100000, 72900000, 76800000,
    77200000,
  ],
  end: [
    3400000, 6500000, 10800000, 16100000, 22700000, 25100000, 27400000, 33500000, 39800000, 40200000, 42800000,
    46800000, 49300000, 52100000, 59500000, 60200000, 63100000, 64600000, 66200000, 69100000, 72900000, 76800000,
    77200000, 83257441,
  ],
  id: [
    'p13.3',
    'p13.2',
    'p13.1',
    'p12',
    'p11.2',
    'p11.1',
    'q11.1',
    'q11.2',
    'q12',
    'q21.1',
    'q21.2',
    'q21.31',
    'q21.32',
    'q21.33',
    'q22',
    'q23.1',
    'q23.2',
    'q23.3',
    'q24.1',
    'q24.2',
    'q24.3',
    'q25.1',
    'q25.2',
    'q25.3',
  ],
  stain: [
    'gneg',
    'gpos50',
    'gneg',
    'gpos75',
    'gneg',
    'acen',
    'acen',
    'gneg',
    'gpos50',
    'gneg',
    'gpos25',
    'gneg',
    'gpos25',
    'gneg',
    'gpos75',
    'gneg',
    'gpos75',
    'gneg',
    'gpos50',
    'gneg',
    'gpos75',
    'gneg',
    'gpos25',
    'gneg',
  ],
};

/** Where the synthetic region starts on chr17: TP53's first position (7,668,402) less the gene's offset in it. */
export const RegionOnChr17 = 7_666_402;

/** A spec with nodes added: data under `data`, navigators under the view's `navigators` (by id, params). */
function withNavigators(
  spec: Node,
  view: string,
  navigators: Readonly<Record<string, NavigatorParams>>,
  data: readonly Node[] = [],
): Node {
  const id = view.slice(view.lastIndexOf('/') + 1);
  const children = (spec.children ?? []).map((c) => {
    if (c.kind === 'data') return { ...c, children: [...(c.children ?? []), ...data] };
    if (c.kind !== 'views') return c;
    return {
      ...c,
      children: (c.children ?? []).map((v) => {
        if (v.id !== id) return v;
        const added = Object.entries(navigators).map(([n, params]) => ({
          kind: 'navigator',
          id: n,
          params: params as unknown as Readonly<Record<string, Value>>,
        }));
        // into the view's navigators when it declares some already (the Globins' overview), else a container of them
        const kids = v.children ?? [];
        const has = kids.some((k) => k.kind === 'navigators');
        return {
          ...v,
          children: has
            ? kids.map((k) => (k.kind === 'navigators' ? { ...k, children: [...(k.children ?? []), ...added] } : k))
            : [...kids, { kind: 'navigators', id: 'navigators', children: added }],
        };
      }),
    };
  });
  return { ...spec, children };
}

/**
 * p53's gene region with chr17's ideogram: the chromosome as a system of its own (its length only), the region placed
 * on it by an offset map, the bands on the chromosome; the region's window within it, the rest washed.
 */
export function geneRegionIdeogram(): Node {
  return withNavigators(
    geneByRegion.build(),
    '/views/main',
    { chr17: ideogramNavigator('/views/main/region', '/data/chr17-bands', { system: '/data/chr17' }) },
    chr17Data(),
  );
}

/**
 * The chromosome as a system of its own (its length and bands from a `seq` source, no letters), the region placed on it
 * by an offset map, and its bands: built (normalized) by the builder, then added to the region's spec.
 */
function chr17Data(): readonly Node[] {
  const b = createSeqSpec();
  const source = chr17BandsSource(b);
  const chr17 = b.data.sequence('chr17', { alphabet: 'dna', source, name: 'chr17' });
  b.data.map('region-on-chr17', { from: RegionSequence, to: chr17, method: 'offset', by: RegionOnChr17 });
  b.data.ranges('chr17-bands', { target: chr17, source });
  return b.build().children?.find((c) => c.kind === 'data')?.children ?? [];
}

/** chr17's bands as an inline `seq` source (`chr17-source`, its system `chr17`): for any spec showing chr17. */
export function chr17BandsSource(b: ReturnType<typeof createSeqSpec>) {
  return b.data.source('chr17-source', { format: 'seq', value: seqToValue(chr17Seq()) as never });
}

/** The level chr17's bands are also served at: a tile of it (1,024 bins) covers the chromosome in two. */
export const Chr17BandsLevel = 65_536;

/**
 * chr17 as a `seq` source: its length (no letters) and its bands, at level 1 and at a coarse level (a navigator reads
 * a ranges layer at its coarse level: the whole chromosome in a tile or two).
 */
export function chr17Seq(): SeqData {
  const bands = {
    form: 'list' as const,
    start: Int32Array.from(Chr17Bands.start),
    end: Int32Array.from(Chr17Bands.end),
    id: categoryColumn(Chr17Bands.id),
    fields: [{ name: 'stain', column: categoryColumn(Chr17Bands.stain) }],
  };
  const pair = (level: number) => ({
    request: { system: 'chr17', start: 0, end: Chr17Length, level, shapes: ['ranges' as const] },
    answer: { ok: true as const, tables: { ranges: bands } },
  });
  return {
    coverage: {
      systems: [{ name: 'chr17', length: Chr17Length, alphabet: 'dna' }],
      shapes: [{ system: 'chr17', shape: 'ranges', fields: [{ name: 'stain', type: 'category' }] }],
      levels: [1, Chr17BandsLevel],
    },
    pairs: [pair(1), pair(Chr17BandsLevel)],
  };
}

/** p53's tracks with a strip of them small above. */
export function tracksStrip(): Node {
  return withNavigators(tracks.build(), '/views/main', {
    strip: stripNavigator('/views/main/protein', { content: 'tracks' }),
  });
}

/** The protein on its gene with a strip of the protein showing its window and its gene follower's. */
export function proteinGeneStrip(): Node {
  return withNavigators(proteinOnGene.build(), '/views/main', {
    strip: stripNavigator(['/views/main/protein', '/views/main/gene']),
  });
}

/** The Globins' 10,000 rows with their vertical overview at the right (the scale view's) and a miniature above. */
export function globinsNavigators(): Node {
  return withNavigators(globins.scale(), '/views/main', {
    miniature: miniatureNavigator('/views/main/rows'),
  });
}

/**
 * The plasmid's ring and linear detail with a strip of the detail's whole system under it (next to what it navigates):
 * the strands' features small, its window wrapping the origin.
 */
export function plasmidStrip(): Node {
  const spec = plasmidRing.build();
  const seq = (spec.children ?? []).find((c) => c.kind === 'data')?.children?.find((c) => c.kind === 'sequence')?.id;
  return withNavigators(spec, '/views/main', {
    strip: stripNavigator('/views/main/detail', {
      system: `/data/${seq}`,
      content: 'tracks',
      tracks: ['/views/main/detail/plus', '/views/main/detail/minus'],
      place: 'bottom',
      size: 16,
    }),
  });
}
