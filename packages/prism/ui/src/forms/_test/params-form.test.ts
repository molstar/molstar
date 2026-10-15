/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsFormModel, paramSections } from '../params-form.js';

const params = {
  type: PD.Select('cartoon', [
    ['cartoon', 'Cartoon'],
    ['surface', 'Surface'],
  ]),
  alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, PD.Essential),
  probe: PD.Numeric(1.4, { min: 0, max: 3 }, { hideIf: (v: any) => v.type !== 'surface' }),
  hidden: PD.Text('x', { isHidden: true }),
};

describe('ParamsFormModel', () => {
  it('holds values, streams each field only when it changes, and reports changes', () => {
    const changes: string[] = [];
    const form = new ParamsFormModel(params, undefined, { onChange: (_v, name) => changes.push(name) });
    const alphas: number[] = [];
    form.field('alpha').subscribe((a) => alphas.push(a));
    form.set('type', 'surface');
    form.set('alpha', 0.5);
    form.set('alpha', 0.5);
    expect(alphas).toEqual([1, 0.5]);
    expect(changes).toEqual(['type', 'alpha']);
    expect(form.values.value).toMatchObject({ type: 'surface', alpha: 0.5 });
  });

  it('follows hideIf on the current values', () => {
    const form = new ParamsFormModel(params);
    const seen: boolean[] = [];
    form.visible('probe').subscribe((v) => seen.push(v));
    form.set('type', 'surface');
    form.set('type', 'cartoon');
    expect(seen).toEqual([false, true, false]);
  });

  it('splits essentials from advanced, skipping hidden params', () => {
    expect(new ParamsFormModel(params).layout()).toEqual({ essentials: ['alpha'], advanced: ['type', 'probe'] });
    const plain = new ParamsFormModel({ a: PD.Boolean(true), b: PD.Text('') });
    expect(plain.layout()).toEqual({ essentials: ['a', 'b'], advanced: [] });
  });

  it('sections params as v5: essentials and the rest, uncategorized first, then categories by name', () => {
    const sections = paramSections({
      a: PD.Numeric(1, {}, { category: 'Shading' }),
      b: PD.Numeric(1, {}, PD.Essential),
      c: PD.Numeric(1, {}, { category: 'Culling' }),
      d: PD.Numeric(1),
      e: PD.Value(1),
      f: PD.Numeric(1, {}, { category: 'Shading' }),
    });
    expect(sections).toEqual({
      essentials: [{ names: ['b'] }],
      advanced: [{ names: ['d'] }, { category: 'Culling', names: ['c'] }, { category: 'Shading', names: ['a', 'f'] }],
    });
  });

  it('keeps a fold state per path for the form', () => {
    const form = new ParamsFormModel(params);
    expect(form.fold('/x#open', true).value).toBe(true);
    form.fold('#advanced').set(true);
    expect(form.fold('#advanced').value).toBe(true);
    expect(new ParamsFormModel(params).fold('#advanced').value).toBe(false);
  });

  it('loads authoritative values silently and resets to defaults', () => {
    const changes: string[] = [];
    const form = new ParamsFormModel(params, undefined, { onChange: (_v, name) => changes.push(name) });
    form.load({ ...form.values.value, alpha: 0.2 });
    expect(changes).toEqual([]);
    form.reset();
    expect(form.values.value.alpha).toBe(1);
    expect(changes).toEqual(['alpha']);
  });
});
