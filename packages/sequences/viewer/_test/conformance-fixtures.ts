/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Documents for the conformance kit's tests (XR1b): the text layout of any document, a matrix, a map row, and sites
// with a K_d each and `custom["org.example"]` for `org.example:stars` and the field provider.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { allTypes } from './prepare-fixtures.js';

/** Every view of a document in the text layout (wrapped). */
export function asText(doc: SeqViewSpecDoc): SeqViewSpecDoc {
  const d = structuredClone(doc);
  for (const v of d.root.children ?? []) if (v.kind === 'view') v.params = { ...v.params, layout: 'text' };
  return d;
}
/** Every section of a document on rings (`style: "circular"`, CG1b). */
export function ringed(doc: SeqViewSpecDoc): SeqViewSpecDoc {
  const d = structuredClone(doc);
  type N = { kind: string; params?: Record<string, unknown>; children?: N[] };
  const walk = (n: N) => {
    if (n.kind === 'section') n.params = { ...n.params, style: 'circular' };
    n.children?.forEach(walk);
  };
  walk(d.root as N);
  return d;
}
/**
 * Every section of a document laid out in `direction` (AX1); `ranged`: each block of a sequence given by value shows
 * its [2, length − 2) only.
 */
export function directed(doc: SeqViewSpecDoc, direction: 1 | -1, ranged = false): SeqViewSpecDoc {
  const d = structuredClone(doc);
  type N = { kind: string; ref?: string; params?: Record<string, unknown>; children?: N[] };
  const lengths = new Map<string, number>();
  const walk = (n: N, f: (n: N) => void) => {
    f(n);
    for (const c of n.children ?? []) walk(c, f);
  };
  walk(d.root as N, (n) => {
    if (n.kind === 'sequence' && n.ref && typeof n.params?.value === 'string')
      lengths.set(n.ref, n.params.value.length);
  });
  walk(d.root as N, (n) => {
    if (n.kind === 'section') n.params = { ...n.params, direction };
    const len = lengths.get(n.params?.sequence as string);
    if (ranged && n.kind === 'block' && len !== undefined && len > 4) n.params = { ...n.params, range: [2, len - 2] };
  });
  return d;
}

/** Every `letters` and `rows` layer of a document drawn through a letter map (LM1): the kit's fixture for it. */
export function mapped(doc: SeqViewSpecDoc, map: string | Record<string, string> = 'complement'): SeqViewSpecDoc {
  const d = structuredClone(doc);
  type N = { kind: string; params?: Record<string, unknown>; children?: N[] };
  const walk = (n: N) => {
    if (n.kind === 'layer' && (n.params?.type === 'letters' || n.params?.type === 'rows'))
      n.params = { ...n.params, letter_map: map };
    n.children?.forEach(walk);
  };
  walk(d.root as N);
  return d;
}

export function matrixDoc() {
  const b = createSeqViewSpec({ title: 'Matrix' });
  const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR'.repeat(10) });
  const index = Array.from({ length: 100 }, (_, i) => i);
  const value = { index, value: index.map((i) => [i % 3, (i * 7) % 5, 1]) };
  const m = b.data.annotation({ ref: 'm', geometry: 'matrix', target: a, states: ['x', 'y', 'z'], value } as never);
  b.view({ ref: 'v', label: 'V' }).block({ sequence: a }).track({ label: 'M' }).matrix({ data: m });
  return b.getDocument();
}
export function mapDoc() {
  const b = createSeqViewSpec({ title: 'Map' });
  const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR'.repeat(10) });
  const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  v.track({ label: 'Map' }).layer({ type: 'seqv:map', fallback: 'omit' });
  return b.getDocument();
}

/**
 * An alignment given by value (8 rows of 90 columns, gaps here and there), its logo and rows: `rows` and `logo` in the
 * kit's edge geometries, which rewrite only what is given by value (case 3 reads its alignment from a file).
 */
export function alignedDoc() {
  const b = createSeqViewSpec({ title: 'Aligned' });
  const base = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYT';
  const rows = Array.from({ length: 8 }, (_, r) =>
    Array.from(base.slice(0, 90), (c, i) => ((i * 7 + r * 13) % 23 === 0 ? '-' : (i + r) % 11 === 0 ? 'A' : c)).join(
      '',
    ),
  );
  const labels = rows.map((_, r) => `r${r}`);
  const aln = b.data.alignment({ ref: 'aln', value: rows, labels, alphabet: 'protein' });
  const profile = b.data.computed({ ref: 'profile', provider: 'profile', alignment: aln });
  const v = b.view({ ref: 'v', label: 'V' }).block({ alignment: aln });
  v.track({ label: 'Logo' }).logo({ data: profile });
  v.track({ label: 'Members' }).rows();
  return b.getDocument();
}

/**
 * Sites with a K_d each and `custom["org.example"]`: stars sized and colored by their affinity, and ranges; `long`:
 * 3,960 positions, the sites spread over them, so the fit is below a px a position (LC-18).
 */
export function starsDoc(o: { ranges?: boolean; rules?: boolean; long?: boolean } = {}) {
  const b = createSeqViewSpec({ title: 'Stars' });
  const a = b.data.sequence({
    ref: 'a',
    label: 'A',
    value: 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'.repeat(o.long ? 120 : 2),
  });
  const start = o.long ? [3, 1200, 2500, 3900] : [3, 12, 30, 47];
  const sites = b.data.annotation({
    ref: 'sites',
    geometry: 'ranges',
    target: a,
    value: {
      start,
      stop: start.map((s) => s + 1),
      kd_nm: [1, 10, 100, 1000],
      label: ['s1', 's2', 's3', 's4'],
    },
    custom: { 'org.example': { unit: 'pKd' } },
  } as never);
  const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  const palette = { kind: 'continuous', colors: ['#000000', '#ffffff'] };
  const color = { theme: 'field_palette', params: { field: 'org.example:affinity', palette } };
  const stars = { type: 'org.example:stars', fallback: 'points', data: sites, color, options: { star_points: 5 } };
  v.track({ label: 'Stars' }).layer(stars as never);
  if (o.ranges) v.track({ label: 'Ranges' }).ranges({ data: sites, color } as never);
  const strong = { when: [{ measure: { field: 'org.example:affinity' }, min: 8 }], color: '#ff0000', label: 'strong' };
  const rules = {
    theme: 'rules',
    params: { rules: [strong], default: { theme: 'uniform', params: { color: '#808080' } } },
  };
  if (o.rules) v.track({ label: 'Rules' }).points({ data: sites, color: rules } as never);
  return b.getDocument();
}

/** A sequence with ranges and values tracks and a `seqv:map` of them (or of the named tracks). */
export function miniatureDoc(tracks?: string[]) {
  const b = createSeqViewSpec({ title: 'Map' });
  const a = b.data.sequence({ ref: 'a', label: 'A', value: 'MKTAYIAKQR'.repeat(40) });
  const hot = b.data.annotation({
    ref: 'hot',
    geometry: 'ranges',
    target: a,
    value: { start: [10, 120, 300], stop: [60, 180, 320], category: ['a', 'b', 'a'] },
  } as never);
  const vals = b.data.annotation({
    ref: 'vals',
    geometry: 'values',
    target: a,
    value: {
      index: Array.from({ length: 400 }, (_, i) => i),
      value: Array.from({ length: 400 }, (_, i) => (i % 40) / 40),
    },
  } as never);
  const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  v.track({ ref: 't.hot', label: 'Hot' }).ranges({ data: hot } as never);
  v.track({ ref: 't.vals', label: 'Vals' }).values({ data: vals } as never);
  v.track({ label: 'Map' }).layer({
    type: 'seqv:map',
    ...(tracks ? { options: { tracks } } : {}),
    fallback: 'omit',
  } as never);
  return b.getDocument();
}

/**
 * `allTypes` with items where edge geometries split them (the kit's `edges`): ranges, points and variants over the
 * cyclic origin's shown seam ([25, 35) at origin 30), a long deletion, a boundary at each end of chain A, and an item
 * across a filtered stretch.
 */
export function edgeTypes(): SeqViewSpecDoc {
  const d = allTypes();
  type N = { kind: string; ref?: string; params?: Record<string, unknown>; children?: N[] };
  const root = d.root as N;
  const data = root.children?.find((c) => c.kind === 'data') as N;
  data.children?.push({
    kind: 'annotation',
    ref: 'edge',
    params: {
      geometry: 'ranges',
      target: 'a',
      value: { start: [25, 100, 0, 300, 40, 61], stop: [35, 290, 0, 300, 41, 66], alt: ['', '', 'M', 'W', '', ''] },
    },
  });
  const section = (root.children?.find((c) => c.kind === 'view')?.children ?? []).find(
    (c) => c.kind === 'section',
  ) as N;
  for (const type of ['ranges', 'points', 'variants', 'emphasis'])
    section.children?.push({
      kind: 'track',
      ref: `edge.${type}`,
      params: { label: `Edge ${type}` },
      children: [{ kind: 'layer', ref: `edge.${type}.layer`, params: { type, data: 'edge' } }],
    });
  return d;
}

/** A document with a `frame` layer over `ra`'s items (labeled), spanning its own track: the kit draws its titles. */
export function withFrame(d: SeqViewSpecDoc): SeqViewSpecDoc {
  type N = { kind: string; ref?: string; params?: Record<string, unknown>; children?: N[] };
  const root = d.root as N;
  const section = (root.children?.find((c) => c.kind === 'view')?.children ?? []).find(
    (c) => c.kind === 'section',
  ) as N;
  section.children?.push({
    kind: 'track',
    ref: 'frame.track',
    params: { label: 'Frames' },
    children: [
      { kind: 'layer', ref: 'frame.layer', params: { type: 'frame', data: 'ra', span: { from: 'frame.track' } } },
    ],
  });
  return d;
}
