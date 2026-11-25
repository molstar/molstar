/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The theme slots in the language (_spec/README.md "Themes", Forms): every form valid as the description reads it, a
// misuse for each `seq.lang.theme.*` code, a string color still valid, and the builder's theme forms (a field over
// inline data listing its values, a ranges layer colored by its category column, a shared scale).

import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import type { Node, Value } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { checkSeqSpec, createSeqSpec, specBounds, theme } from '../index.js';

function spec(color: unknown, more: { scale?: boolean; type?: string; channel?: string } = {}) {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEPQS' });
  const r = b.data.ranges('r', {
    target: s,
    rows: { start: new Int32Array([0, 2]), end: new Int32Array([2, 5]), category: ['a', 'b'] },
  });
  const v = b.data.values('v', { target: s, values: new Float32Array([1, 2, 3, 4, 5, 6]) });
  if (more.scale) b.views.scale('cons', { domain: [0, null] as unknown as [number, number] });
  b.views.view('main', (view) => {
    const t = view.section('p', { of: s }).track('t');
    const type = more.type ?? 'ranges';
    t.layer('l', {
      type: type as 'ranges',
      ...(type === 'letters' ? {} : { data: type === 'values' ? v : r }),
      [more.channel ?? 'color']: color as Value,
    });
  });
  return b.build();
}

const codes = (n: Node, bounds = specBounds()) => checkSeqSpec(n, bounds).map((i) => i.code);

describe('the theme slots', () => {
  it('a string color, a token, light-dark and every form are valid', () => {
    for (const c of [
      '#d62728',
      'token:range',
      'light-dark(#000000, token:range)',
      theme.uniform('#ff0000'),
      theme.field('category'),
      theme.scheme('categories', { field: 'category' }),
      theme.rules([theme.rule([theme.eq('category', 'a')], '#ff0000')], 'token:missing'),
    ])
      expect(codes(spec(c)), JSON.stringify(c)).toEqual([]);
    expect(codes(spec(theme.field('score', { scale: '/views/cons' }), { scale: true }))).toEqual([]);
    expect(
      codes(
        spec(theme.letters({ colors: { C: '#e0c000' }, groups: [theme.group('AILMV', '#ffafaf', 'aliphatic')] }), {
          type: 'letters',
        }),
      ),
    ).toEqual([]);
    expect(codes(spec(0.5, { channel: 'size' }))).toEqual([]);
    expect(codes(spec('dashed', { channel: 'line' }))).toEqual([]);
  });

  it('prints a slot in place as an S-expression', () => {
    const s = printSexpr(spec(theme.field('score', { domain: [0, 1] })));
    expect(s).toContain('(field :name score :domain [0 1])');
  });

  it('a misuse for each code', () => {
    expect(codes(spec('reddish'))).toEqual(['seq.lang.theme.color']);
    expect(codes(spec(3))).toEqual(['seq.lang.theme.channel']);
    expect(codes(spec('#ff0000', { channel: 'size' }))).toEqual(['seq.lang.theme.channel']);
    expect(codes(spec('wiggly', { channel: 'line' }))).toEqual(['seq.lang.theme.channel']);
    expect(codes(spec(theme.continuous({ colors: 'viridis' })))).toEqual(['seq.lang.theme.form']);
    expect(
      codes(
        spec(theme.field('start', { palette: theme.binned({ edges: [0, 10, 5], colors: ['#000000', '#ffffff'] }) })),
      ),
    ).toEqual(['seq.lang.theme.stops']);
    expect(codes(spec(theme.field('start', { domain: [5, 1] })))).toEqual(['seq.lang.theme.stops']);
    expect(codes(spec(theme.field('score', { scale: '/views/nope' })))).toContain('seq.lang.theme.scale');
    expect(codes(spec(theme.rules([theme.rule([theme.eq('category', 3)], '#ff0000')])))).toEqual([
      'seq.lang.theme.type',
    ]);
    expect(codes(spec(theme.rules([theme.rule([theme.lt('category', 3)], '#ff0000')])))).toEqual([
      'seq.lang.theme.type',
    ]);
    let deep: Node = theme.uniform('#ff0000');
    for (let i = 0; i < 6; i++) deep = theme.rules([], deep);
    expect(codes(spec(deep))).toEqual(['seq.lang.theme.depth']);
    expect(codes(spec(deep), specBounds({ themeDepth: 8 }))).toEqual([]);
  });

  it('a non-string non-form value is refused by the description', () => {
    expect(codes(spec({ nope: 1 })).some((c) => c.startsWith('repr.'))).toBe(true);
  });
});

describe('the builder’s theme forms', () => {
  it('a ranges layer over rows with a category is colored by it, its values listed in first-occurrence order', () => {
    const n = spec(undefined);
    const layer = JSON.stringify(n);
    expect(layer).toContain('"palette":{"kind":"categorical","params":{"colors":"categories","values":["a","b"]}}');
  });

  it('theme.field over inline text lists the values, unless asked not to; numbers are never listed', () => {
    expect(JSON.stringify(spec(theme.field('category', { list: false })))).not.toContain('"values":["a","b"]');
    expect(JSON.stringify(spec(theme.field('start')))).not.toContain('"values":["0"');
    expect(JSON.stringify(spec(theme.field('category')))).not.toContain('__list');
  });

  it('b.views.scale adds a scale node under views', () => {
    const n = spec('#ff0000', { scale: true });
    const views = n.children!.find((c) => c.kind === 'views')!;
    expect(views.children!.map((c) => `${c.kind}:${c.id}`)).toContain('scale:cons');
  });
});
