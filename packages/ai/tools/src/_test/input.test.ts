/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ObjectShape } from '../index.js';
import { shapeInput } from '../index.js';

const props = (shape: ObjectShape) =>
  (shapeInput(shape).jsonSchema as { properties: Record<string, unknown> }).properties;
const message = (shape: ObjectShape, args: Record<string, unknown>) => {
  const r = shapeInput(shape).parse(args);
  return r.ok ? undefined : r.message;
};

describe('a field shape, as the JSON Schema the model sees and the check (tools-and-session.md §2)', () => {
  it('is an object schema in draft 2020-12, required where not optional and without a default', () => {
    const shape = {
      kind: 'object',
      fields: {
        scene: { kind: 'string', description: 'Scene id.' },
        note: { kind: 'string', optional: true },
        size: { kind: 'number', default: 2 },
      },
    } as const;
    expect(shapeInput(shape).jsonSchema).toEqual({
      $schema: 'https://json-schema.org/draft/2020-12/schema',
      type: 'object',
      properties: {
        scene: { type: 'string', description: 'Scene id.' },
        note: { type: 'string' },
        size: { type: 'number', default: 2 },
      },
      required: ['scene'],
    });
    expect(shapeInput(shape).parse({ scene: 's' })).toEqual({ ok: true, value: { scene: 's', size: 2 } });
  });

  it('string: an enum, a format said in the description', () => {
    const shape = {
      kind: 'object',
      fields: {
        repr: { kind: 'string', enum: ['cartoon', 'ball-and-stick'], description: 'Representation.' },
        target: { kind: 'string', format: 'A:315, B:200-210', description: 'Where.' },
        free: { kind: 'string', format: 'a name' },
      },
    } as const;
    expect(props(shape)).toEqual({
      repr: { type: 'string', enum: ['cartoon', 'ball-and-stick'], description: 'Representation.' },
      target: { type: 'string', description: 'Where. Format: A:315, B:200-210.' },
      free: { type: 'string', description: 'Format: a name.' },
    });
    expect(message(shape, { repr: 'x', target: 't', free: 'f' })).toBe(
      'repr: expected one of cartoon, ball-and-stick (got "x").',
    );
    expect(message(shape, { target: 't', free: 'f' })).toBe(
      'repr: required, one of cartoon, ball-and-stick (Representation.).',
    );
    expect(shapeInput(shape).parse({ repr: 'cartoon', target: 't', free: 'f' }).ok).toBe(true);
  });

  it('number and integer: bounds', () => {
    const shape = {
      kind: 'object',
      fields: {
        size: { kind: 'number', min: 0.5, max: 4, description: 'Size factor.' },
        count: { kind: 'integer', min: 1 },
      },
    } as const;
    expect(props(shape)).toEqual({
      size: { type: 'number', minimum: 0.5, maximum: 4, description: 'Size factor.' },
      count: { type: 'integer', minimum: 1 },
    });
    expect(shapeInput(shape).parse({ size: 1, count: 2 }).ok).toBe(true);
    expect(message(shape, { size: 9, count: 2 })).toBe('size: Too big: expected number to be <=4 (Size factor.).');
    expect(message(shape, { size: 1, count: 1.5 })).toBe('count: expected int (got 1.5).');
    expect(message(shape, { size: 'big', count: 2 })).toBe('size: expected number (got "big") (Size factor.).');
    expect(message(shape, { count: 2 })).toBe('size: required (Size factor.).');
  });

  it('boolean, and array of an item with a size limit', () => {
    const shape = {
      kind: 'object',
      fields: {
        label: { kind: 'boolean' },
        frame: { kind: 'array', item: { kind: 'string' }, max: 2, description: 'Targets.' },
      },
    } as const;
    expect(props(shape)).toEqual({
      label: { type: 'boolean' },
      frame: { type: 'array', items: { type: 'string' }, maxItems: 2, description: 'Targets.' },
    });
    expect(message(shape, { label: 'yes', frame: [] })).toBe('label: expected boolean (got "yes").');
    expect(message(shape, { label: true, frame: ['a', 1] })).toBe('frame[1]: expected string (got 1).');
    expect(message(shape, { label: true, frame: ['a', 'b', 'c'] })).toContain('frame: Too big');
  });

  it('object: nested fields, the whole path in a message', () => {
    const shape = {
      kind: 'object',
      fields: {
        claims: {
          kind: 'array',
          item: { kind: 'object', fields: { text: { kind: 'string' }, page: { kind: 'integer' } } },
        },
        box: { kind: 'object', optional: true, fields: { x: { kind: 'number' } } },
      },
    } as const;
    expect(props(shape).claims).toEqual({
      type: 'array',
      items: {
        type: 'object',
        properties: { text: { type: 'string' }, page: { type: 'integer' } },
        required: ['text', 'page'],
      },
    });
    expect(message(shape, { claims: [{ text: 'x', page: 1 }, { text: 'y' }] })).toBe('claims[1].page: required.');
    expect(message(shape, { claims: [], box: { x: 'a' } })).toBe('box.x: expected number (got "a").');
  });

  it('variant: a discriminated union on its tag', () => {
    const shape = {
      kind: 'object',
      fields: {
        mark: {
          kind: 'variant',
          tag: 'type',
          description: 'What to draw.',
          variants: {
            box: { kind: 'object', fields: { size: { kind: 'number' } } },
            arrow: { kind: 'object', fields: { from: { kind: 'string' }, to: { kind: 'string', optional: true } } },
          },
        },
      },
    } as const;
    const mark = props(shape).mark as { oneOf: { properties: Record<string, unknown>; required: string[] }[] };
    expect(mark.oneOf.map((o) => o.properties.type)).toEqual([
      { type: 'string', const: 'box' },
      { type: 'string', const: 'arrow' },
    ]);
    expect(mark.oneOf[1]?.required).toEqual(['type', 'from']);
    expect(shapeInput(shape).parse({ mark: { type: 'arrow', from: 'A' } })).toEqual({
      ok: true,
      value: { mark: { type: 'arrow', from: 'A' } },
    });
    expect(message(shape, { mark: { type: 'box' } })).toBe('mark.size: required.');
    expect(message(shape, { mark: { type: 'cone' } })).toBe('mark.type: expected one of box, arrow (got "cone").');
    expect(message(shape, {})).toBe('mark: required (What to draw.).');
  });

  it('json: any value, required unless optional', () => {
    const shape = {
      kind: 'object',
      fields: { value: { kind: 'json', description: 'A value.' }, more: { kind: 'json', optional: true } },
    } as const;
    expect(props(shape)).toEqual({ value: { description: 'A value.' }, more: {} });
    for (const value of [1, 'a', null, { a: [1] }, [true]]) expect(shapeInput(shape).parse({ value }).ok).toBe(true);
    expect(message(shape, {})).toBe('value: required (A value.).');
  });

  it('a variant without variants is refused when the tool is defined', () => {
    expect(() => shapeInput({ kind: 'object', fields: { m: { kind: 'variant', tag: 't', variants: {} } } })).toThrow(
      'has no variants',
    );
  });
});
