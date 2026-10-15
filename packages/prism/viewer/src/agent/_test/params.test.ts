/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { normalize, shapeInput } from '@molstar/ai-tools';
import { Color } from '@molstar/core/util/color/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { colorOf, paramsShape } from '../params.js';

const Params = {
  repr: PD.Select('cartoon', [
    ['cartoon', 'Cartoon'],
    ['ball-and-stick', 'Ball & Stick'],
  ]),
  size: PD.Numeric(2, { min: 1, max: 10, step: 1 }),
  alpha: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
  label: PD.Text('', { label: 'Label', description: 'What the component is called.' }),
  visible: PD.Boolean(true),
  color: PD.Color(Color(0xff0000)),
  axes: PD.MultiSelect(
    ['x'],
    [
      ['x', 'X'],
      ['y', 'Y'],
    ],
  ),
  range: PD.Interval([0, 1]),
  center: PD.Vec3([0, 0, 0] as any),
  style: PD.Group({
    opacity: PD.Numeric(1, { min: 0, max: 1 }),
    quality: PD.Select('auto', [
      ['auto', 'Auto'],
      ['high', 'High'],
    ]),
  }),
  theme: PD.MappedStatic('uniform', {
    uniform: PD.Group({ value: PD.Color(Color(0xffffff)) }),
    'chain-id': PD.Group({
      palette: PD.Select('a', [
        ['a', 'A'],
        ['b', 'B'],
      ]),
    }),
    none: PD.EmptyGroup(),
  }),
  items: PD.ObjectList({ name: PD.Text(''), weight: PD.Numeric(1) }, (e) => e.name),
  query: PD.Script({ language: 'pymol', expression: 'all' }),
  file: PD.File(),
  secret: PD.Text('x', { isHidden: true }),
};

describe('paramsShape', () => {
  const ps = paramsShape(Params);

  it('maps each kind to a field shape, every field optional', () => {
    const f = ps.shape.fields;
    expect(f.repr).toMatchObject({
      kind: 'string',
      enum: ['cartoon', 'ball-and-stick'],
      optional: true,
      default: 'cartoon',
    });
    expect(f.size).toMatchObject({ kind: 'integer', min: 1, max: 10, default: 2 });
    expect(f.alpha).toMatchObject({ kind: 'number', min: 0, max: 1 });
    expect(f.label).toMatchObject({ kind: 'string', description: 'What the component is called.' });
    expect(f.visible).toMatchObject({ kind: 'boolean', default: true });
    expect(f.color).toMatchObject({ kind: 'string', format: 'color', default: 'red' });
    expect(f.axes).toMatchObject({ kind: 'array', item: { kind: 'string', enum: ['x', 'y'] } });
    expect(f.range).toMatchObject({ kind: 'array', max: 2 });
    expect(f.center).toMatchObject({ kind: 'array', max: 3 });
    expect(f.style).toMatchObject({
      kind: 'object',
      fields: { opacity: { kind: 'number' }, quality: { kind: 'string' } },
    });
    expect(f.theme).toMatchObject({ kind: 'variant', tag: 'name' });
    expect(Object.keys((f.theme as any).variants)).toEqual(['uniform', 'chain-id', 'none']);
    expect(f.items).toMatchObject({ kind: 'array', item: { kind: 'object' } });
    expect(f.query).toMatchObject({ kind: 'object', fields: { expression: { kind: 'string', verbatim: true } } });
    for (const field of Object.values(f)) expect(field.optional).toBe(true);
  });

  it('leaves out what JSON cannot carry, and hidden params', () => {
    expect(ps.omitted).toEqual(['file', 'secret']);
    expect(ps.shape.fields.file).toBeUndefined();
  });

  it('names select options by their labels where they say more', () => {
    expect(ps.shape.fields.repr.description).toBe('Options: ball-and-stick (Ball & Stick).');
  });

  it('makes a JSON Schema and a check through ai-tools', () => {
    const input = shapeInput(ps.shape);
    expect((input.jsonSchema as any).properties.repr.enum).toEqual(['cartoon', 'ball-and-stick']);
    expect(input.parse(normalize(input.shape, { repr: 'BALL-AND-STICK', size: '3' }))).toEqual({
      ok: true,
      value: expect.objectContaining({ repr: 'ball-and-stick', size: 3 }),
    });
    expect(input.parse({ size: 20 }).ok).toBe(false);
    expect(input.parse({ theme: { name: 'chain-id', params: { palette: 'b' } } }).ok).toBe(true);
  });

  it('turns arguments into values over the defaults', () => {
    const values = ps.toValues({
      repr: 'ball-and-stick',
      color: '#00ff00',
      style: { opacity: 0.3 },
      theme: { name: 'chain-id', params: { palette: 'b' } },
      items: [{ name: 'a' }],
      query: { expression: 'resn HEM' },
    });
    expect(values.repr).toBe('ball-and-stick');
    expect(values.size).toBe(2);
    expect(values.color).toBe(Color(0x00ff00));
    expect(values.style).toEqual({ opacity: 0.3, quality: 'auto' });
    expect(values.theme).toEqual({ name: 'chain-id', params: { palette: 'b' } });
    expect(values.items).toEqual([{ name: 'a', weight: 1 }]);
    expect(values.query).toEqual({ language: 'pymol', expression: 'resn HEM' });
    expect(values.file).toBe(null);
    expect(values.secret).toBe('x');
  });

  it("a mapped choice without params takes that choice's defaults", () => {
    expect(ps.toValues({ theme: { name: 'uniform' } }).theme).toEqual({
      name: 'uniform',
      params: { value: Color(0xffffff) },
    });
  });

  it('reads values back as arguments', () => {
    const args = ps.toArgs(PD.getDefaultValues(Params));
    expect(args.color).toBe('red');
    expect(args.theme).toEqual({ name: 'uniform', params: { value: 'white' } });
    expect(args.query).toEqual({ language: 'pymol', expression: 'all' });
    expect(args.file).toBeUndefined();
    expect(ps.toValues(args)).toEqual({ ...PD.getDefaultValues(Params) });
  });
});

describe('colorOf', () => {
  it('reads names, hex in its forms and numbers', () => {
    expect(colorOf('red')).toBe(Color(0xff0000));
    expect(colorOf('Light Blue')).toBe(Color(0xadd8e6));
    expect(colorOf('#0f0')).toBe(Color(0x00ff00));
    expect(colorOf('0x0000ff')).toBe(Color(0x0000ff));
    expect(colorOf('00ff00')).toBe(Color(0x00ff00));
    expect(colorOf(255)).toBe(Color(255));
    expect(colorOf('not a color')).toBeUndefined();
  });
});
