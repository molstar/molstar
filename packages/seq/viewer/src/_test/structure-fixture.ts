/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A structure shaped as 2OCJ for the structure link's tests: P04637 (synthetic letters) and the entity it maps to by a
// `structure` map (residues 94 to 312), four chains A–D of it modeled 96–289 (the ends unobserved), a zinc each, and
// optionally a second map of the same pair one residue off (`aln`, ranked lower). A view of P04637 (its letters and a
// domain at 144–153) above a section of the complex; and a recording side, the structure viewer as its host
// translates it, with its own events to fire.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { location } from '@molstar/seqx-core/query';
import type { StructureSide, StructureTarget } from '../link-structure.js';

export const Up = '/data/p53';
export const Ent = '/data/ent';
export const Complex = '/data/2ocj';
export const Sifts = '/data/sifts';
export const Aln = '/data/aln';
export const UpSection = '/views/main/p53';
export const ComplexSection = '/views/main/complex';
export const UpLength = 393;
/** The entity: P04637 94 to 312 (0-based 93). */
export const EntStart = 93;
export const EntLength = 219;
/** Each chain's modeled residues on the entity. */
export const Observed: readonly [number, number] = [2, 196];
export const Chains = ['A', 'B', 'C', 'D'] as const;

const Letters = 'ACDEFGHIKLMNPQRSTVWY';
const letters = (n: number) => Array.from({ length: n }, (_, i) => Letters[(i * 7 + 3) % 20]).join('');

export function p53Structure(o: { readonly aln?: boolean } = {}): Node {
  const b = createSeqSpec({ title: '2OCJ-shaped' });
  const up = b.data.sequence('p53', { alphabet: 'protein', letters: letters(UpLength) });
  const ent = b.data.sequence('ent', {
    alphabet: 'protein',
    letters: letters(UpLength).slice(EntStart, EntStart + EntLength),
  });
  const zn = b.data.molecule('zn', { class: 'ion', residues: ['ZN'] });
  b.data.complex('2ocj', { name: '2OCJ' }, (c) => {
    for (const chain of Chains)
      c.unit(chain, { of: ent, chain, auth_chain: chain, observed: { start: [Observed[0]], end: [Observed[1]] } });
    for (const chain of Chains) c.unit(`Z${chain}`, { of: zn, chain: `Z${chain}`, auth_chain: chain });
  });
  b.data.map('sifts', {
    from: up,
    to: ent,
    method: 'structure',
    segments: { from: [EntStart], to: [0], length: [EntLength], provenance: ['sifts'] },
    ...(o.aln ? { priority: 1 } : {}),
  });
  if (o.aln)
    b.data.map('aln', {
      from: up,
      to: ent,
      method: 'segments',
      segments: { from: [EntStart + 1], to: [0], length: [EntLength - 1] },
    });
  const domains = b.data.ranges('domains', {
    target: up,
    rows: { start: [EntStart + 50], end: [EntStart + 60], id: ['D'] },
  });
  b.views.view('main', (v) => {
    const s = v.section('p53', { of: up });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('domains').layer('domains', { type: 'ranges', data: domains });
    v.section('complex', { of: Complex }, (s) => s.track('letters').layer('letters', { type: 'letters' }));
  });
  return b.build();
}

/** Positions on P04637 as a Location. */
export const onUp = (start: number, end: number) => location({ parts: [{ system: Up, ranges: [start, end] }] });

/** Positions on each named chain as a target (the side's words for a residue range). */
export const onChains = (start: number, end: number, chains: readonly string[] = Chains) =>
  location({ parts: chains.map((c) => ({ system: `${Complex}/${c}`, ranges: [start, end] })) });

/** A message the side was told, in order. */
export type Told =
  | { readonly kind: 'highlight' | 'select'; readonly target: StructureTarget | null }
  | { readonly kind: 'focus'; readonly target: StructureTarget; readonly durationMs: number };

/** A structure viewer's side that records what it is told and fires its own events by hand. */
export interface RecordingSide extends StructureSide {
  readonly told: Told[];
  readonly abort: AbortController;
  hover(target: StructureTarget | null, at?: number): void;
  click(target: StructureTarget | null): void;
  doubleClick(target: StructureTarget): void;
  /** How many listeners the link left subscribed. */
  readonly listeners: number;
}

export function recordingSide(): RecordingSide {
  const told: Told[] = [];
  const abort = new AbortController();
  const hovers = new Set<(t: StructureTarget | null, at?: number) => void>();
  const selects = new Set<(t: StructureTarget | null) => void>();
  const focuses = new Set<(t: StructureTarget) => void>();
  const sub = <F>(set: Set<F>, fn: F) => {
    set.add(fn);
    return () => void set.delete(fn);
  };
  return {
    told,
    abort,
    signal: abort.signal,
    highlight: (target) => void told.push({ kind: 'highlight', target }),
    select: (target) => void told.push({ kind: 'select', target }),
    focus: (target, { durationMs }) => void told.push({ kind: 'focus', target, durationMs }),
    onHover: (fn) => sub(hovers, fn),
    onSelect: (fn) => sub(selects, fn),
    onFocus: (fn) => sub(focuses, fn),
    hover: (t, at) => {
      for (const f of hovers) f(t, at);
    },
    click: (t) => {
      for (const f of selects) f(t);
    },
    doubleClick: (t) => {
      for (const f of focuses) f(t);
    },
    get listeners() {
      return hovers.size + selects.size + focuses.size;
    },
  };
}
