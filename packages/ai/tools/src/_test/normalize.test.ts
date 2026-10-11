/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { normalize } from '../index.js';

describe('normalizing what small models send, against the shape (tools-and-session.md §2)', () => {
  const shape = {
    kind: 'object',
    fields: {
      scene: { kind: 'string' },
      repr: { kind: 'string', enum: ['cartoon', 'ball-and-stick'] },
      size: { kind: 'number' },
      count: { kind: 'integer' },
      label: { kind: 'boolean' },
      frame: { kind: 'array', item: { kind: 'string' } },
      sizes: { kind: 'array', item: { kind: 'number' } },
      kinds: { kind: 'array', item: { kind: 'string', enum: ['a', 'b'] } },
      residues: { kind: 'string' },
      code: { kind: 'string', verbatim: true },
      replace: { kind: 'string', verbatim: true },
      claims: {
        kind: 'array',
        item: { kind: 'object', fields: { text: { kind: 'string' }, page: { kind: 'integer' } } },
      },
      box: { kind: 'object', fields: { x: { kind: 'number' }, note: { kind: 'string' } } },
      mark: {
        kind: 'variant',
        tag: 'type',
        variants: {
          box: { kind: 'object', fields: { size: { kind: 'number' } } },
          label: { kind: 'object', fields: { on: { kind: 'boolean' }, text: { kind: 'string' } } },
        },
      },
      any: { kind: 'json' },
    },
  } as const;
  const n = (args: Record<string, unknown>) => normalize(shape, args);

  it('an empty string or null is absent', () => {
    expect(n({ scene: '', size: null, note: '  ', label: 'true' })).toEqual({ label: true });
  });

  it('strings are trimmed, but not verbatim ones, and an empty verbatim string is a value', () => {
    expect(n({ scene: ' pocket ', code: '  x();\n', replace: '' })).toEqual({
      scene: 'pocket',
      code: '  x();\n',
      replace: '',
    });
  });

  it('"3" is a number and "true" a boolean where one is asked for; a number is a string where a string is', () => {
    expect(n({ size: ' 2.5 ', count: '3', label: 'false', residues: 315, scene: 7 })).toEqual({
      size: 2.5,
      count: 3,
      label: false,
      residues: '315',
      scene: '7',
    });
    // What does not read as asked is left for the check to refuse.
    expect(n({ size: 'big', label: 'yes' })).toEqual({ size: 'big', label: 'yes' });
  });

  it('a list of strings is one comma-joined text where one text is asked for', () => {
    expect(n({ residues: ['A:315', ' A:318'] })).toEqual({ residues: 'A:315, A:318' });
  });

  it('"A, B" is a list where a list of non-objects is', () => {
    expect(n({ frame: 'A:315, site.pocket,', sizes: '1, 2.5', kinds: 'A, b' })).toEqual({
      frame: ['A:315', 'site.pocket'],
      sizes: [1, 2.5],
      kinds: ['a', 'b'],
    });
    expect(n({ frame: [' a ', 'b'] })).toEqual({ frame: ['a', 'b'] });
  });

  it('an enum matches in any case', () => {
    expect(n({ repr: 'Ball-And-Stick' })).toEqual({ repr: 'ball-and-stick' });
    expect(n({ repr: 'spheres' })).toEqual({ repr: 'spheres' });
  });

  it('objects in a list are read by their own fields, and so are objects', () => {
    expect(
      n({ claims: [{ text: ' found ', page: '3', gone: '' }, 'odd'], box: { x: '1', note: ' hi ', y: ' ' } }),
    ).toEqual({
      claims: [{ text: 'found', page: 3 }, 'odd'],
      box: { x: 1, note: 'hi' },
    });
    // A string is not split into objects.
    expect(n({ claims: 'a, b' })).toEqual({ claims: 'a, b' });
  });

  it('a variant is read by the fields of the variant its tag names', () => {
    expect(n({ mark: { type: ' Box ', size: '2' } })).toEqual({ mark: { type: 'box', size: 2 } });
    expect(n({ mark: { type: 'label', on: 'true', text: 5 } })).toEqual({
      mark: { type: 'label', on: true, text: '5' },
    });
    // No such variant: left as sent, for the check to name the choices.
    expect(n({ mark: { type: 'cone', size: '2' } })).toEqual({ mark: { type: 'cone', size: '2' } });
  });

  it('a json value is kept, trimmed if it is text', () => {
    expect(n({ any: ' x ' })).toEqual({ any: 'x' });
    expect(n({ any: { a: ' ' } })).toEqual({ any: { a: ' ' } });
  });

  it('keeps fields the shape does not name, trimmed', () => {
    expect(n({ other: ' x ' })).toEqual({ other: 'x' });
  });
});
