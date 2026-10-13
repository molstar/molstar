/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { interpolate, interpolateIdPath, validateIdPathTemplate } from '../string.js';

describe('string', () => {
  describe('interpolate', () => {
    it('substitutes simple placeholders', () => {
      expect(interpolate('Hello ${name}!', { name: 'world' })).toBe('Hello world!');
    });

    it('substitutes multiple placeholders', () => {
      expect(interpolate('${a}/${b}', { a: 'x', b: 'y' })).toBe('x/y');
    });

    it('uses empty string for missing keys', () => {
      expect(interpolate('${missing}', {})).toBe('');
    });

    it('does not evaluate expressions', () => {
      expect(interpolate('${triggers}', { triggers: '<i>click</i>' })).toBe('<i>click</i>');
    });
  });

  describe('interpolateIdPath', () => {
    it('substitutes id', () => {
      expect(interpolateIdPath('./data/${id}.bcif', '1abc')).toBe('./data/1abc.bcif');
    });

    it('substitutes id.toLowerCase()', () => {
      expect(interpolateIdPath('./${id.toLowerCase()}.mdb', 'EMD-1234')).toBe('./emd-1234.mdb');
    });

    it('substitutes id.toUpperCase()', () => {
      expect(interpolateIdPath('./${id.toUpperCase()}.mdb', 'emd-1234')).toBe('./EMD-1234.mdb');
    });

    it('substitutes id.substr()', () => {
      expect(interpolateIdPath('./${id.substr(1, 2)}/${id}.bcif', '1abc')).toBe('./ab/1abc.bcif');
    });

    it('substitutes id.substring() with an exclusive end index', () => {
      expect(interpolateIdPath('./${id.substring(1, 2)}/${id}.bcif', '1abc')).toBe('./a/1abc.bcif');
      expect(interpolateIdPath('./${id.substring(1, 3)}/${id}.bcif', '1abc')).toBe('./ab/1abc.bcif');
    });

    it('swaps reversed substring bounds', () => {
      expect(interpolateIdPath('${id.substring(3, 1)}', '1abc')).toBe('ab');
    });

    it('handles empty and out-of-range substrings', () => {
      expect(interpolateIdPath('${id.substring(2, 2)}', '1abc')).toBe('');
      expect(interpolateIdPath('${id.substring(2, 20)}', '1abc')).toBe('bc');
      expect(interpolateIdPath('${id.substr(20, 2)}', '1abc')).toBe('');
      expect(interpolateIdPath('${id.substr(1, 0)}', '1abc')).toBe('');
    });

    it('clamps length to available characters', () => {
      expect(interpolateIdPath('./${id.substr(0, 20)}.bcif', '1234567890')).toBe('./1234567890.bcif');
      expect(interpolateIdPath('./${id.substr(8, 20)}.bcif', '1234567890')).toBe('./90.bcif');
    });

    it('substitutes id.substring() with one arg', () => {
      expect(interpolateIdPath('./${id.substring(1)}.bcif', 'x1abc')).toBe('./1abc.bcif');
    });

    it('substitutes id.slice() with one argument', () => {
      expect(interpolateIdPath('./${id.slice(1)}.bcif', 'x1abc')).toBe('./1abc.bcif');
    });

    it('substitutes id.slice() with two arguments', () => {
      expect(interpolateIdPath('./${id.slice(1, 3)}.bcif', '1abc')).toBe('./ab.bcif');
    });

    it.each(['substr(1, 2)', 'substring(1, 3)', 'substring(1)', 'slice(1, 3)', 'slice(1)'])(
      'applies case conversion before and after %s',
      (call) => {
        const expected = call.endsWith('(1)') ? 'abc' : 'ab';
        expect(interpolateIdPath('${id.toLowerCase().' + call + '}', '1ABC')).toBe(expected);
        expect(interpolateIdPath('${id.' + call + '.toLowerCase()}', '1ABC')).toBe(expected);
        expect(interpolateIdPath('${id.toUpperCase().' + call + '}', '1abc')).toBe(expected.toUpperCase());
        expect(interpolateIdPath('${id.' + call + '.toUpperCase()}', '1abc')).toBe(expected.toUpperCase());
      },
    );

    it('applies multiple calls to the result of the previous call', () => {
      expect(
        interpolateIdPath('${id.toUpperCase().substring(1, 4).substr(1, 2).slice(0, 1).toLowerCase()}', '1abcd'),
      ).toBe('b');
    });

    it('preserves call order when case conversion changes string length', () => {
      expect(interpolateIdPath('${id.toUpperCase().slice(0, 2)}', 'ßabc')).toBe('SS');
      expect(interpolateIdPath('${id.slice(0, 2).toUpperCase()}', 'ßabc')).toBe('SSA');
    });

    it('preserves the original id across placeholders', () => {
      expect(interpolateIdPath('./${id.substring(1, 3).toLowerCase()}/${id}.bcif', '1ABC')).toBe('./ab/1ABC.bcif');
    });

    it('handles config template patterns', () => {
      const template = './path-to-binary-cif/${id.substr(1, 2)}/${id}.bcif';
      expect(interpolateIdPath(template, '1abc')).toBe('./path-to-binary-cif/ab/1abc.bcif');
    });

    it('rejects unsupported expressions', () => {
      expect(() => interpolateIdPath('./${id + 1}.bcif', '1abc')).toThrow('Unsupported id path expression');
    });
  });

  describe('validateIdPathTemplate', () => {
    it('accepts supported templates', () => {
      expect(() => validateIdPathTemplate('./${id.substr(1, 2)}/${id}.bcif')).not.toThrow();
      expect(() =>
        validateIdPathTemplate('./${id.substring(1, 3).toLowerCase()}/${id.toLowerCase()}.bcif'),
      ).not.toThrow();
      expect(() => validateIdPathTemplate('./${id.toUpperCase().slice(1, 3).toLowerCase()}.bcif')).not.toThrow();
    });

    it('rejects unsupported templates', () => {
      expect(() => validateIdPathTemplate('./${id + 1}.bcif')).toThrow('Unsupported id path expression');
    });
  });

  it.each([
    'other.toLowerCase()',
    'idOther.toLowerCase()',
    'id.toLowerCase(1)',
    'id.toUpperCase(1)',
    'id.substr(1)',
    'id.substring()',
    'id.slice(1, 2, 3)',
    'id.slice(-1)',
    'id.slice(1.5)',
    'id.slice(1 + 1)',
    'id.toLowerCase().replace("a", "b")',
    'id.toLowerCase().constructor()',
    'id.toLowerCase().slice(1) + id',
    'id.toLowerCase()trailing',
    'id.toLowerCase().',
  ])('rejects unsupported expressions in both validation and interpolation: %s', (expression) => {
    const template = '${' + expression + '}';
    expect(() => validateIdPathTemplate(template)).toThrow('Unsupported id path expression');
    expect(() => interpolateIdPath(template, '1ABC')).toThrow('Unsupported id path expression');
  });
});
