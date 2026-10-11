/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { argsFor, formFields, formOnlyParams, keptValues, readValue, valueLabel } from '../flat-fields.js';

const params = {
  label: PD.Text('', { label: 'Label' }),
  representation: PD.Select('cartoon', [
    ['cartoon', 'Cartoon'],
    ['spacefill', 'Spacefill'],
  ]),
  wholeResidues: PD.Boolean(true),
  radius: PD.Numeric(5, { min: 0, max: 20, step: 0.5 }),
  count: PD.Numeric(3),
  level: PD.Select(1, [
    [1, 'One'],
    [2, 'Two'],
  ]),
  color: PD.Color(0xff0000 as any),
  options: PD.Group({ quality: PD.Text('auto') }),
  secret: PD.Text('x', { isHidden: true }),
  note: PD.Optional(PD.Text('', { label: 'Note' })) as PD.Any,
};

describe('formFields', () => {
  const fields = formFields(params, { representation: { before: 'as' }, label: { format: 'name' } });

  it('reads simple params as inline fields in order, with labels, hints and what is required', () => {
    expect(fields.map((f) => [f.id, f.kind, f.label, f.required])).toEqual([
      ['label', 'text', 'Label', true],
      ['representation', 'select', 'Representation', true],
      ['wholeResidues', 'boolean', 'Whole Residues', true],
      ['radius', 'slider', 'Radius', true],
      ['count', 'number', 'Count', true],
      ['note', 'text', 'Note', false],
    ]);
    expect(fields[0]?.format).toBe('name');
    expect(fields[1]?.before).toBe('as');
    expect(fields[3]).toMatchObject({ min: 0, max: 20, step: 0.5 });
  });

  it('leaves non-string selects, colors and groups to the full form, and hidden params out', () => {
    expect(formOnlyParams(params, fields)).toEqual(['level', 'color', 'options', 'secret']);
  });
});

describe('values', () => {
  const fields = formFields(params);
  const field = (id: string) => fields.find((f) => f.id === id)!;

  it('reads values as their fields hold them, or rejects them', () => {
    expect(readValue(field('representation'), 'spacefill')).toBe('spacefill');
    expect(readValue(field('representation'), 'sticks')).toBeUndefined();
    expect(readValue(field('wholeResidues'), 'false')).toBe(false);
    expect(readValue(field('wholeResidues'), true)).toBe(true);
    expect(readValue(field('wholeResidues'), 'maybe')).toBeUndefined();
    expect(readValue(field('radius'), 21)).toBeUndefined();
    expect(readValue(field('radius'), Number.NaN)).toBeUndefined();
    expect(readValue(field('count'), 1e6)).toBe(1e6);
    expect(readValue(field('label'), '')).toBeUndefined();
    expect(readValue(field('note'), '')).toBe('');
  });

  it('keeps only values that fit', () => {
    expect(keptValues(fields, { representation: 'sticks', radius: 2, wholeResidues: 'true', unknown: 'x' })).toEqual({
      radius: 2,
      wholeResidues: true,
    });
  });

  it('shows values by their option label', () => {
    expect(valueLabel(field('representation'), 'spacefill')).toBe('Spacefill');
    expect(valueLabel(field('wholeResidues'), false)).toBe('off');
    expect(valueLabel(field('radius'), 2.5)).toBe('2.5');
  });

  it("gives a run the inline values over the params' defaults; optional params left empty stay absent", () => {
    const args = argsFor(params, { label: 'Pocket', radius: 2 });
    expect(args).toMatchObject({ label: 'Pocket', radius: 2, representation: 'cartoon', wholeResidues: true });
    expect(args.options).toEqual({ quality: 'auto' });
    expect('note' in args).toBe(false);
    expect(argsFor(params, { note: 'hi' }).note).toBe('hi');
  });
});
