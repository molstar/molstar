/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { interpolate, interpolateIdPath, validateIdPathTemplate } from '../string';

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

        it('substitutes id.substring() with two args using substr-like length semantics', () => {
            expect(interpolateIdPath('./${id.substring(1, 2)}/${id}.bcif', '1abc')).toBe('./ab/1abc.bcif');
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
        });

        it('rejects unsupported templates', () => {
            expect(() => validateIdPathTemplate('./${id + 1}.bcif')).toThrow('Unsupported id path expression');
        });
    });
});
