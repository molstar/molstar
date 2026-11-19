/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Rings and blocks fixtures: a small cyclic construct (600 nt, a feature over its origin, one on the minus strand) as
// a ring with a linear detail following its window, and a small complex's chains and ligands as blocks.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { seeded } from './fixtures.js';

export const RingLength = 600;

/** The construct's letters: seeded, so drawing is checked against them. */
export function ringLetters(): string {
  const r = seeded(11);
  let s = '';
  for (let i = 0; i < RingLength; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s;
}

/** Its features: `ori` over the origin (start past end), `amp` on the minus strand, `mcs` short. */
export const RingFeatures = {
  start: [560, 120, 300],
  end: [40, 420, 330],
  id: ['ori', 'amp', 'mcs'],
  label: ['ori', 'AmpR', 'MCS'],
  strand: ['+', '-', '+'],
};

export interface RingOptions {
  readonly rotation?: number;
  /** The linear detail follows the ring (true by default); with `ringFollows`, the ring follows the detail. */
  readonly detail?: boolean;
  readonly ringFollows?: boolean;
  readonly labels?: 'outside' | 'none';
}

/** The construct as a ring (features, letters, a ruler) and its linear detail. */
export function ringSpec(options: RingOptions = {}): Node {
  const b = createSeqSpec({ title: 'ring' });
  const c = b.data.sequence('construct', { alphabet: 'dna', letters: ringLetters(), cyclic: true });
  const f = b.data.ranges('features', { target: c, rows: RingFeatures });
  b.views.view('v', (v) => {
    const ringParams = {
      of: c,
      layout: 'rings' as const,
      ...(options.rotation !== undefined ? { rotation: options.rotation } : {}),
      ...(options.labels ? { labels: options.labels } : {}),
    };
    const detailFirst = options.ringFollows === true;
    const detail = detailFirst ? v.section('detail', { of: c }) : undefined;
    const ring = v.section('ring', detail ? { ...ringParams, follow: detail } : ringParams);
    ring.track('features').layer('features', { type: 'ranges', data: f });
    ring.track('ruler').layer('ruler', { type: 'ruler' });
    if (detail) {
      detail.track('features').layer('features', { type: 'ranges', data: f });
      detail.track('letters').layer('letters', { type: 'letters' });
    } else if (options.detail !== false) {
      const d = v.section('detail', { of: c, follow: ring });
      d.track('features').layer('features', { type: 'ranges', data: f });
      d.track('letters').layer('letters', { type: 'letters' });
    }
  });
  return b.build();
}

/**
 * A plasmid-sized ring (2,686 nt): four rings (features on each strand, sites, a ruler), about 30 labels, a linear
 * detail following its window; what the ring's frame target is measured on.
 */
export function plasmidRingSpec(): Node {
  const n = 2686;
  const r = seeded(23);
  let letters = '';
  for (let i = 0; i < n; i++) letters += 'ACGT'[Math.floor(r() * 4)];
  const plus = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
  const minus = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
  for (let k = 0; k < 13; k++) {
    const s = Math.floor((k * n) / 13 + r() * 40);
    const f = k % 2 ? minus : plus;
    f.start.push(s);
    f.end.push((s + 80 + Math.floor(r() * 120)) % n);
    f.id.push(`f${k}`);
    f.label.push(`feature ${k}`);
  }
  // one across the origin on each strand
  plus.start.push(2600);
  plus.end.push(60);
  plus.id.push('ori');
  plus.label.push('ori');
  minus.start.push(2500);
  minus.end.push(120);
  minus.id.push('bla');
  minus.label.push('bla');
  const sites = { start: [] as number[], end: [] as number[], id: [] as string[], label: [] as string[] };
  for (let k = 0; k < 14; k++) {
    const s = Math.floor(r() * (n - 6));
    sites.start.push(s);
    sites.end.push(s + 6);
    sites.id.push(`s${k}`);
    sites.label.push(`Site${k}`);
  }
  const b = createSeqSpec({ title: 'plasmid ring' });
  const c = b.data.sequence('plasmid', { alphabet: 'dna', letters, cyclic: true });
  const fp = b.data.ranges('plus', { target: c, rows: plus });
  const fm = b.data.ranges('minus', { target: c, rows: minus });
  const fs = b.data.ranges('sites', { target: c, rows: sites });
  b.views.view('v', (v) => {
    const ring = v.section('ring', { of: c, layout: 'rings' });
    ring.track('plus').layer('plus', { type: 'ranges', data: fp });
    ring.track('minus').layer('minus', { type: 'ranges', data: fm });
    ring.track('sites', { height: 12 }).layer('sites', { type: 'ranges', data: fs });
    ring.track('ruler').layer('ruler', { type: 'ruler' });
    const d = v.section('detail', { of: c, follow: ring });
    d.track('plus').layer('plus', { type: 'ranges', data: fp });
    d.track('minus').layer('minus', { type: 'ranges', data: fm });
    d.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** A small complex: two chains (one in two copies), a ligand, as blocks of one section. */
export function blocksSpec(options: { arrange?: 'row' | 'stack'; copies?: 'merge' | 'each'; layout?: 'rings' } = {}) {
  const b = createSeqSpec({ title: 'blocks' });
  const r = seeded(5);
  const protein = (n: number) => {
    let s = '';
    for (let i = 0; i < n; i++) s += 'ACDEFGHIKLMNPQRSTVWY'[Math.floor(r() * 20)];
    return s;
  };
  const a = b.data.sequence('cdk2', { alphabet: 'protein', letters: protein(120) });
  const c = b.data.sequence('cyclin', { alphabet: 'protein', letters: protein(80) });
  const l = b.data.sequence('atp', { alphabet: 'protein', letters: 'X' });
  const fa = b.data.ranges('cdk2-features', {
    target: a,
    rows: { start: [40], end: [52], id: ['psta'], label: ['PSTAIRE'] },
  });
  const fc = b.data.ranges('cyclin-features', {
    target: c,
    rows: { start: [10], end: [60], id: ['box'], label: ['cyclin box'] },
  });
  b.views.view('v', (v) => {
    const s = v.section('complex', {
      ...(options.arrange ? { arrange: options.arrange } : {}),
      ...(options.copies ? { copies: options.copies } : {}),
      ...(options.layout ? { layout: options.layout } : {}),
    });
    s.block('a', { of: a, label: 'CDK2' });
    s.block('b', { of: c, label: 'cyclin A' });
    s.block('c', { of: a, label: 'CDK2' });
    s.block('atp', { of: l, label: 'ATP' });
    s.track('features', {}, (t) => {
      t.layer('cdk2', { type: 'ranges', data: fa });
      t.layer('cyclin', { type: 'ranges', data: fc });
    });
    s.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}
