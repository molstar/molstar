/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The documents of `Seq*/Color schemes` (lib/seqstar/storybook/src/_stories/color-schemes.stories.tsx), kept out of the story
// file so the smoke test (lib/seqstar/storybook/_test/smoke.test.ts) can open each on a headless viewer.

import { createSeqViewSpec, type SchemeRef } from '@molstar/seq-core/builder';
import { type ColorSchemeEntry, KernelColorSchemes, type SeqViewSpecDoc } from '@molstar/seq-core/format';

type B = ReturnType<typeof createSeqViewSpec>;
const label = (name: string) => KernelColorSchemes.find((e) => e.name === name)?.label ?? name;
const names = (kind: (e: ColorSchemeEntry) => boolean) => KernelColorSchemes.filter(kind).map((e) => e.name);
export const PROTEIN = names((e) => 'alphabets' in e.scheme && e.scheme.alphabets.includes('protein'));
export const NUCLEIC = names((e) => 'alphabets' in e.scheme && e.scheme.alphabets.includes('rna'));
const VALUES = names((e) => e.scheme.kind === 'sequential' || e.scheme.kind === 'diverging');
const CATEGORIES = names((e) => e.scheme.kind === 'categorical');

/** One section per scheme over the same alignment, each following the first. */
export function alignmentGallery(title: string, file: string, schemes: readonly string[]): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title });
  const aln = b.data.alignment({ ref: 'aln', source: { host_file: file, format: 'stockholm' } });
  const v = b.view({ ref: 'v', label: title });
  let first: ReturnType<ReturnType<B['view']>['section']> | null = null;
  for (const name of schemes) {
    const sec = v.section({ ref: `v.${name}`, label: label(name), ...(first ? { follow: first.handle } : {}) });
    sec
      .block({ alignment: aln })
      .track({ label: label(name) })
      .rows({ ref: `rows.${name}`, color: b.colors.scheme(name) });
    first ??= sec;
  }
  return b.getDocument();
}

const SEQ =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNT';
const at = <T>(n: number, f: (i: number) => T) => Array.from({ length: n }, (_, i) => f(i));
/** p53 with a score in 0..1 (sequential), a change in -1..1 (diverging) and seven categories of regions. */
export function valuesGallery(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Values and categories' });
  const p = b.data.sequence({ ref: 'p53', label: 'p53', value: SEQ });
  const n = SEQ.length;
  const index = at(n, (i) => i);
  const score = b.data.annotation({
    ref: 'score',
    geometry: 'values',
    target: p,
    value: { index, value: at(n, (i) => Number((0.5 + 0.5 * Math.sin(i / 23)).toFixed(3))) },
  });
  const change = b.data.annotation({
    ref: 'change',
    geometry: 'values',
    target: p,
    value: { index, value: at(n, (i) => Number(Math.sin(i / 31).toFixed(3))) },
  });
  const kinds = ['helix', 'strand', 'loop', 'binding', 'disordered', 'motif', 'domain'];
  const regions = b.data.annotation({
    ref: 'regions',
    geometry: 'ranges',
    target: p,
    value: {
      start: at(22, (k) => k * 14),
      stop: at(22, (k) => k * 14 + 12),
      category: at(22, (k) => kinds[k % 7]),
    },
  });
  const v = b.view({ ref: 'v', label: 'Values and categories' });
  const sec = v.section({ ref: 'v.main', label: 'p53' }).block({ sequence: p });
  sec.track({ label: '' }).ruler();
  for (const name of VALUES) {
    const data = name === 'rd_bu' ? change : score;
    sec.track({ label: label(name) }).color({ ref: `values.${name}`, data, color: b.colors.scheme(name) });
  }
  for (const name of CATEGORIES)
    sec.track({ label: label(name) }).ranges({ ref: `cats.${name}`, data: regions, color: b.colors.scheme(name) });
  return b.getDocument();
}

/** A host's own scheme, registered: listed by the picker after the built-ins. */
export const PASTEL = {
  name: 'story:pastel',
  label: "Pastel groups (a host's)",
  description: "a host's own scheme, registered at runtime",
  provenance: 'this story',
  colorblind: false,
  scheme: {
    kind: 'letters',
    alphabets: ['protein'],
    groups: [
      { letters: 'AVLIMFWC', color: '#cfe3f7', label: 'hydrophobic' },
      { letters: 'DE', color: '#f7cfd8', label: 'acidic' },
      { letters: 'KRH', color: '#d4cff7', label: 'basic' },
      { letters: 'STNQGPY', color: '#d8f0d0', label: 'other' },
    ],
  },
} as const;
export const HOST = [{ svs: { colorSchemes: [PASTEL] } }];

export function pickerDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Pick a scheme' });
  const aln = b.data.alignment({ ref: 'sh2', source: { host_file: 'PF00017.seed.sto', format: 'stockholm' } });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
  const v = b.view({ ref: 'v', label: 'SH2' }).block({ alignment: aln });
  v.track({ label: 'Conservation' }).color({ ref: 'cons.color', data: cons });
  v.track({ label: 'Members' }).rows({ ref: 'members' });
  return b.getDocument();
}

/**
 * 25 tracks of three kinds on p53, as track templates give them: five genes with nothing to list, ten
 * transcripts by part (exon, UTR) and ten domains by kind. Each kind's legends are identical.
 */
export function manyDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Many tracks' });
  const p = b.data.sequence({ ref: 'p53', label: 'p53', value: SEQ });
  const ranges = (ref: string, start: number[], len: number, category?: string[]) =>
    b.data.annotation({
      ref,
      geometry: 'ranges',
      target: p,
      value: { start, stop: start.map((x) => x + len), ...(category ? { category } : {}) },
    });
  const sec = b
    .view({ ref: 'v', label: 'Many tracks' })
    .section({ ref: 'v.main', label: 'p53' })
    .block({ sequence: p });
  sec.track({ label: '' }).ruler();
  for (let g = 1; g <= 5; g++) {
    sec.track({ label: `Gene G${g}` }).ranges({ data: ranges(`gene${g}`, [g * 40], 90) });
    for (let k = 1; k <= 2; k++) {
      const T = `G${g}.${k}`;
      const parts = ranges(`tx${T}`, [g * 40, g * 40 + 30, g * 40 + 60], 20, ['UTR', 'exon', 'exon']);
      sec.track({ label: `Transcript ${T}` }).ranges({ data: parts });
      const dom = ranges(`dom${T}`, [g * 40 + 10 * k], 25, ['domain']);
      sec.track({ label: `${T} domain` }).ranges({ data: dom });
    }
  }
  return b.getDocument();
}

const custom = (b: B): Record<string, SchemeRef> => ({
  'custom letters: groups and a default': b.colors.custom({
    kind: 'letters',
    alphabets: ['protein'],
    groups: [
      { letters: 'DE', color: '#d62728', label: 'acidic' },
      { letters: 'KRH', color: '#1f77b4', label: 'basic' },
      { letters: 'C', color: '#e6b800', label: 'cysteine' },
    ],
    default: '#d9d9d9',
    dark: { default: '#3a3d44' },
  }),
  'derived: Zappo with cysteine magenta': b.colors.from(b.colors.zappo, { overrides: { C: '#ff00cc' } }),
});
export function customDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Your own palette' });
  const aln = b.data.alignment({ ref: 'sh2', source: { host_file: 'PF00017.seed.sto', format: 'stockholm' } });
  const cons = b.data.computed({ ref: 'cons', provider: 'column_stats', alignment: aln, measure: 'conservation' });
  const v = b.view({ ref: 'v', label: 'Custom' });
  const ramp = b.colors.custom({ kind: 'sequential', colors: ['#fff7e6', '#f4a259', '#8c2d04'] });
  let first: ReturnType<ReturnType<B['view']>['section']> | null = null;
  for (const [title, color] of Object.entries(custom(b))) {
    const sec = v.section({ label: title, ...(first ? { follow: first.handle } : {}) }).block({ alignment: aln });
    if (!first) sec.track({ label: 'Custom ramp' }).color({ data: cons, color: ramp });
    sec.track({ label: 'Members' }).rows({ color });
    first ??= sec;
  }
  return b.getDocument();
}

export type Shared = 'own' | 'group' | 'explicit';
/**
 * Two signals on p53 with different ranges, a score in 0..1 and a count in 0..8, as bars (a `values` layer) and as color
 * (a `color` layer over Viridis). `own`: each layer's domain is its data's. `group`: all four name one `domain_group`.
 * `explicit`: the same group, and the score's color gives its own `domain` (0..16).
 */
export function sharedDoc(mode: Shared): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Shared scales' });
  const p = b.data.sequence({ ref: 'p53', label: 'p53', value: SEQ });
  const n = SEQ.length;
  const column = (ref: string, f: (i: number) => number) =>
    b.data.annotation({
      ref,
      geometry: 'values',
      target: p,
      value: { index: at(n, (i) => i), value: at(n, (i) => Number(f(i).toFixed(3))) },
    });
  const score = column('score', (i) => 0.5 + 0.5 * Math.sin(i / 23));
  const count = column('count', (i) => 4 + 4 * Math.sin(i / 37 + 1));
  // A track's `height` is in rows (22 px): bars 1.3 rows, color strips 0.6; the section as tall as its rows. A `values`
  // layer has no color of its own (the ink: black), so the bars take the scheme too.
  const sec = b
    .view({ ref: 'v', label: 'Shared scales' })
    .section({ ref: 'v.main', label: 'p53', height: 'content' })
    .block({ sequence: p });
  sec.track({ label: '' }).ruler();
  const group = mode === 'own' ? {} : { domain_group: 'signal' };
  const viridis = (domain?: [number, number]) => b.colors.viridis(domain ? { domain } : {});
  const scoreColor = viridis(mode === 'explicit' ? [0, 16] : undefined);
  const countColor = viridis();
  const bars = 1.3;
  const strip = 0.6;
  sec
    .track({ label: 'Score bars', height: bars })
    .values({ ref: 'score.bars', data: score, color: scoreColor, ...group });
  sec
    .track({ label: 'Score color', height: strip })
    .color({ ref: 'score.color', data: score, color: scoreColor, ...group });
  sec
    .track({ label: 'Count bars', height: bars })
    .values({ ref: 'count.bars', data: count, color: countColor, ...group });
  sec
    .track({ label: 'Count color', height: strip })
    .color({ ref: 'count.color', data: count, color: countColor, ...group });
  return b.getDocument();
}

/** One signal in 0..8 as color strips, each over its own range: the data's, wider than it, and narrower than it. */
export function rangeDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Color range' });
  const p = b.data.sequence({ ref: 'p53', label: 'p53', value: SEQ });
  const n = SEQ.length;
  const count = b.data.annotation({
    ref: 'count',
    geometry: 'values',
    target: p,
    value: {
      index: at(n, (i) => i),
      value: at(n, (i) => Number((4 + 4 * Math.sin(i / 37 + 1)).toFixed(3))),
    },
  });
  const sec = b
    .view({ ref: 'v', label: 'Color range' })
    .section({ ref: 'v.main', label: 'p53', height: 'content' })
    .block({ sequence: p });
  sec.track({ label: '' }).ruler();
  const strips: [string, string, [number, number] | undefined][] = [
    ['own', "The data's range, 0 to 8", undefined],
    ['wider', 'Wider, 0 to 20', [0, 20]],
    ['narrower', 'Narrower, 3 to 5', [3, 5]],
  ];
  for (const [ref, label, domain] of strips)
    sec
      .track({ label, height: 0.8 })
      .color({ ref: `range.${ref}`, data: count, color: b.colors.viridis(domain ? { domain } : {}) });
  return b.getDocument();
}
