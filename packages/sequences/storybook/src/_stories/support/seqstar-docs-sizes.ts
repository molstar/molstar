/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The documents of `Seq*/Themes/Size themes` (lib/seqstar/storybook/src/_stories/size-themes.stories.tsx), kept out of the story
// file so the smoke test (lib/seqstar/storybook/_test/smoke.test.ts) can open each on a headless viewer.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

const SEQ =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNT';
const at = <T>(n: number, f: (i: number) => T) => Array.from({ length: n }, (_, i) => f(i));
const scale = (field: string, range: [number, number], more: Record<string, unknown> = {}) => ({
  theme: 'field_scale',
  params: { field, range, ...more },
});

/**
 * Four tracks sized by a field of their data: contacts stroked by how often they were seen (width) with the same count
 * as color, ranges as thick as their confidence, points as large as their score, bars as wide as their value.
 */
export function sizesDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Size themes' });
  const p = b.data.sequence({ ref: 'p53', label: 'p53', value: SEQ });
  const counts = [12, 3, 7, 1, 9, 2, 5, 11, 4, 6];
  const from = [5, 30, 55, 75, 100, 120, 150, 170, 200, 230];
  const contacts = b.data.annotation({
    ref: 'contacts',
    geometry: 'pairs',
    target: p,
    value: {
      start: from,
      b_target: from.map(() => 'p53'),
      b_start: from.map((s, i) => s + 20 + 6 * i),
      category: from.map((_, i) => (i % 2 === 0 ? 'plus' : 'minus')),
      count: counts,
    },
  } as never);
  const regions = b.data.annotation({
    ref: 'regions',
    geometry: 'ranges',
    target: p,
    value: {
      start: [20, 70, 130, 200, 260],
      stop: [60, 120, 185, 250, 300],
      confidence: [0.2, 0.9, 0.5, 1, 0.35],
      score: [8, 90, 30, 100, 15],
      category: ['a', 'b', 'c', 'd', 'e'],
    },
  } as never);
  const signal = b.data.annotation({
    ref: 'signal',
    geometry: 'values',
    target: p,
    value: { index: at(50, (i) => i * 6), value: at(50, (i) => Number((5 + 4 * Math.sin(i / 5)).toFixed(2))) },
  });
  const v = b.view({ ref: 'v', label: 'Size themes' });
  const sec = v.section({ ref: 'v.main', label: 'p53', height: 'content' }).block({ sequence: p });
  sec.track({ label: '' }).ruler();
  sec.track({ label: 'Contacts', height: 4 }).pairs({
    data: contacts,
    width: scale('count', [1, 8]),
    radius: scale('count', [2, 5]),
    color: { theme: 'field_palette', params: { field: 'count', palette: { kind: 'continuous', colors: 'Blues' } } },
  } as never);
  sec.track({ label: 'Confidence', height: 1.4 }).ranges({
    data: regions,
    thickness: scale('confidence', [0.2, 1]),
  } as never);
  sec.track({ label: 'Score', height: 1.6 }).points({ data: regions, radius: scale('score', [2, 7]) } as never);
  sec.track({ label: 'Value', height: 2 }).values({
    data: signal,
    color: b.colors.viridis(),
    thickness: scale('value', [0.3, 1]),
  } as never);
  return b.getDocument();
}

/** The same score through each scale into the same range, as the radius of a point: linear, sqrt and log. */
export function scalesDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Scales' });
  const p = b.data.sequence({ ref: 'p53', label: 'p53', value: SEQ });
  const sites = b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: p,
    value: {
      start: at(12, (i) => 10 + i * 25),
      stop: at(12, (i) => 11 + i * 25),
      score: [1, 2, 4, 8, 16, 32, 64, 100, 50, 25, 10, 3],
      category: at(12, () => 'site'),
    },
  } as never);
  const sec = b
    .view({ ref: 'v', label: 'Scales' })
    .section({ ref: 'v.main', label: 'p53', height: 'content' })
    .block({ sequence: p });
  sec.track({ label: '' }).ruler();
  for (const kind of ['linear', 'sqrt', 'log'])
    sec
      .track({ label: kind, height: 1.5 })
      .points({ ref: `sites.${kind}`, data: sites, radius: scale('score', [2, 8], { scale: kind }) } as never);
  return b.getDocument();
}

/**
 * A sashimi plot (`b.preset.sashimi`): a synthetic 900-nucleotide gene with four exons, its RNA-seq coverage as bars,
 * and its splice junctions as arcs stroked and labeled by read count, the + strand's above the line and the - strand's
 * (two antisense junctions) below. An exon-skipping junction (exon 1 to 3) is thin: few reads.
 */
export function sashimiDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'Sashimi' });
  const bases = at(900, (i) => 'ACGT'[(i * 7 + (i >> 3)) % 4]).join('');
  const g = b.data.sequence({ ref: 'gene', label: 'GENE1', alphabet: 'dna', value: bases });
  const exons: [number, number][] = [
    [50, 180],
    [300, 420],
    [520, 640],
    [760, 880],
  ];
  const exon = b.data.annotation({
    ref: 'exons',
    geometry: 'ranges',
    target: g,
    value: {
      start: exons.map(([s]) => s),
      stop: exons.map(([, e]) => e),
      label: exons.map((_, i) => `exon ${i + 1}`),
      category: exons.map(() => 'exon'),
    },
  } as never);
  // Coverage: high over the exons, a little intronic noise, a dip where exon 2 is skipped.
  const depth = at(900, (i) => {
    const k = exons.findIndex(([s, e]) => i >= s && i < e);
    return k < 0 ? (i % 37 === 0 ? 2 : 0) : [44, 39, 35, 30][k] + ((i * 13) % 7) - 3;
  });
  const coverage = b.data.annotation({
    ref: 'coverage',
    geometry: 'values',
    target: g,
    value: { index: at(900, (i) => i), value: depth },
  });
  const junctions = b.data.annotation({
    ref: 'junctions',
    geometry: 'pairs',
    target: g,
    value: {
      start: [179, 419, 179, 639, 700, 400],
      b_target: at(6, () => 'gene'),
      b_start: [300, 520, 520, 760, 560, 250],
      count: [42, 37, 6, 30, 4, 2],
      strand: ['+', '+', '+', '+', '-', '-'],
    },
  } as never);
  const sec = b
    .view({ ref: 'v', label: 'Sashimi' })
    .section({ ref: 'v.main', label: 'GENE1', height: 'content' })
    .block({ sequence: g });
  sec.track({ label: '' }).ruler();
  sec.track({ label: 'Exons' }).ranges({ data: exon, labels: true });
  sec.track({ label: 'Coverage', height: 2 }).values({ data: coverage });
  b.preset.sashimi(sec, junctions, { label: 'Junctions', range: [1, 8] });
  return b.getDocument();
}
