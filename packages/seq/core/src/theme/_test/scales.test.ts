/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Palettes through their scales (_spec/README.md "Themes", Scales): a log scale's non-positive values missing (or
// clamped to the low end) in binned and continuous palettes alike; MolViewSpec's palettes converted with their scale;
// a center making the default palette diverge; the sizes legend on the channel's own scale; the ramp's ticks where the
// scale puts them; only the lists the schemes name.

import { LightTokens } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import {
  builtInSchemes,
  compileChannel,
  compilePalette,
  fromMvsPalette,
  fromMvsScale,
  legendOf,
  namedList,
  type ThemeBody,
  type ThemeNode,
  toMvsPalette,
} from '../index.js';

const schemes = builtInSchemes();
const env = (over: Record<string, unknown> = {}) => ({ tokens: LightTokens, schemes, ...over });
const binned: ThemeNode = { kind: 'binned', params: { bins: 3, colors: 'ramp' } };
const ramp: ThemeNode = { kind: 'continuous', params: { colors: 'ramp' } };

describe('a log scale', () => {
  it('a value at or below zero is missing, in a binned palette too (never the first bin)', () => {
    const p = compilePalette(binned, { tokens: LightTokens, domain: [1, 1000], scale: { type: 'log' } });
    expect(p.num(0)).toBeUndefined();
    expect(p.num(-5)).toBeUndefined();
    expect(p.num(0.5)).toBe(p.num(1));
  });

  it('`non_positive: clamp` draws it at the low end, binned and continuous, and MolViewSpec carries it', () => {
    for (const node of [binned, ramp]) {
      const scale = { type: 'log' as const, non_positive: 'clamp' as const };
      const p = compilePalette(node, { tokens: LightTokens, domain: [1, 1000], scale });
      expect(p.num(0)).toBe(p.num(1));
      expect(p.num(-5)).toBe(p.num(1));
      expect(toMvsPalette(node, LightTokens, { ...scale, domain: [1, 1000] })).toMatchObject({ non_positive: 'clamp' });
    }
  });
});

describe('MolViewSpec’s palettes with their scale', () => {
  it('round-trip a log scale, a domain, a center and the colors past the ends', () => {
    const node: ThemeNode = { kind: 'continuous', params: { colors: 'ramp', under: '#000000', over: '#ffffff' } };
    const scale = { type: 'log' as const, domain: [1, 1000] as const, mid: 10 };
    const mvs = toMvsPalette(node, LightTokens, scale)!;
    expect(mvs).toMatchObject({ scale: 'log', mid: 10, value_domain: [1, 1000] });
    const back = fromMvsPalette(mvs);
    const backScale = fromMvsScale(mvs);
    expect(backScale).toEqual({ type: 'log', domain: [1, 1000], mid: 10 });
    expect(toMvsPalette(back, LightTokens, backScale)).toEqual(mvs);
    const a = compilePalette(node, { tokens: LightTokens, domain: [1, 1000], scale });
    const b = compilePalette(back, { tokens: LightTokens, domain: [1, 1000], scale: backScale });
    for (const v of [0.5, 1, 3, 10, 70, 1000, 2000]) expect(b.num(v), String(v)).toBe(a.num(v));
  });
});

describe('a center', () => {
  it('makes a numeric field’s default palette diverge: the center is the diverging middle', () => {
    const ch = compileChannel(
      'color',
      { kind: 'field', params: { name: 'v', domain: [-1, 1], mid: 0 } } as ThemeBody,
      'token:range',
      env({ numeric: () => true }),
    );
    const at = (v: number) => <number | undefined>ch.at({ value: v, letter: '', start: 0, end: 0, field: () => v });
    // the center within the lookup table's step of the middle token
    const near = (a: number, b: number) =>
      [24, 16, 8].every((k) => Math.abs(((a >>> k) & 255) - ((b >>> k) & 255)) <= 6);
    expect(near(at(0)!, LightTokens.colors['seq-div-mid'])).toBe(true);
    expect(at(-1)).toBe(LightTokens.colors['seq-div-neg']);
    expect(at(1)).toBe(LightTokens.colors['seq-div-pos']);
  });
});

describe('legends on the channel’s scale', () => {
  it('the sizes legend is on the radius’s sqrt scale, as the marks are', () => {
    const ch = compileChannel(
      'size',
      { kind: 'field', params: { name: 'f', domain: [0, 100], range: [0, 1] } } as ThemeBody,
      1,
      env({ numeric: () => true, radius: true }),
    );
    const l = legendOf(ch, { title: 'f', layer: 'l' })!;
    expect(l.entries.kind).toBe('sizes');
    if (l.entries.kind !== 'sizes') return;
    for (const { value, size } of l.entries.items)
      expect(size).toBeCloseTo(ch.at({ value, letter: '', start: 0, end: 0, field: () => value }) as number, 6);
  });

  it('a ramp’s ticks are where the scale puts them: a log scale’s decades, a center at the middle', () => {
    const log = compileChannel(
      'color',
      { kind: 'field', params: { name: 'v', domain: [1, 1000], type: 'log' } } as ThemeBody,
      'token:range',
      env({ numeric: () => true }),
    );
    const l = legendOf(log, { title: 'v', layer: 'l' })!;
    if (l.entries.kind !== 'ramp') throw new Error('a ramp');
    expect(l.entries.ticks).toEqual([1, 10, 100, 1000]);
    expect(l.entries.tickAt[1]).toBeCloseTo(1 / 3, 6);
    const mid = compileChannel(
      'color',
      { kind: 'field', params: { name: 'v', domain: [0, 100], mid: 80 } } as ThemeBody,
      'token:range',
      env({ numeric: () => true }),
    );
    const m = legendOf(mid, { title: 'v', layer: 'l' })!;
    if (m.entries.kind !== 'ramp') throw new Error('a ramp');
    expect(m.entries.tickAt[m.entries.ticks.indexOf(80)]).toBe(0.5);
  });
});

describe('named lists', () => {
  it('only the lists the schemes name, by this language’s names', () => {
    for (const n of [
      'viridis',
      'cividis',
      'magma',
      'blues',
      'greens',
      'greens-r',
      'rd-bu',
      'okabe-ito',
      'tableau10',
      'category10',
    ])
      expect(namedList(n, LightTokens), n).toBeDefined();
    // greens-r is greens, dark to light
    expect(namedList('greens-r', LightTokens)).toEqual([...namedList('greens', LightTokens)!].reverse());
    expect(namedList('plasma', LightTokens)).toBeUndefined();
    expect(namedList('red-blue', LightTokens)).toBeUndefined();
  });
});
