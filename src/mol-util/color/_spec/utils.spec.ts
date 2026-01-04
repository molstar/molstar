/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Zachary Charlop-Powers
 */

import { Color } from '../color';
import { decodeColor, decodeColorWithAlpha } from '../utils';

describe('color utils', () => {
    describe('decodeColor', () => {
        it('6-digit hex codes', () => {
            expect(decodeColor('#ff0000')).toBe(Color.fromRgb(255, 0, 0));
            expect(decodeColor('#00ff00')).toBe(Color.fromRgb(0, 255, 0));
            expect(decodeColor('#0000ff')).toBe(Color.fromRgb(0, 0, 255));
            expect(decodeColor('#ffffff')).toBe(Color.fromRgb(255, 255, 255));
            expect(decodeColor('#000000')).toBe(Color.fromRgb(0, 0, 0));
        });

        it('6-digit hex codes (uppercase)', () => {
            expect(decodeColor('#FF0000')).toBe(Color.fromRgb(255, 0, 0));
            expect(decodeColor('#00FF00')).toBe(Color.fromRgb(0, 255, 0));
            expect(decodeColor('#0000FF')).toBe(Color.fromRgb(0, 0, 255));
        });

        it('3-digit hex codes', () => {
            expect(decodeColor('#f00')).toBe(Color.fromRgb(255, 0, 0));
            expect(decodeColor('#0f0')).toBe(Color.fromRgb(0, 255, 0));
            expect(decodeColor('#00f')).toBe(Color.fromRgb(0, 0, 255));
            expect(decodeColor('#fff')).toBe(Color.fromRgb(255, 255, 255));
        });

        it('8-digit hex codes (strips alpha)', () => {
            // Alpha channel should be stripped, only RGB returned
            expect(decodeColor('#ff0000ff')).toBe(Color.fromRgb(255, 0, 0));
            expect(decodeColor('#00ff0080')).toBe(Color.fromRgb(0, 255, 0));
            expect(decodeColor('#0000ff00')).toBe(Color.fromRgb(0, 0, 255));
            expect(decodeColor('#ffffff7f')).toBe(Color.fromRgb(255, 255, 255));
        });

        it('named colors', () => {
            expect(decodeColor('red')).toBe(0xff0000);
            expect(decodeColor('blue')).toBe(0x0000ff);
            expect(decodeColor('white')).toBe(0xffffff);
        });

        it('rgb() format', () => {
            expect(decodeColor('rgb(255, 0, 0)')).toBe(Color.fromRgb(255, 0, 0));
            expect(decodeColor('rgb(0, 255, 0)')).toBe(Color.fromRgb(0, 255, 0));
            expect(decodeColor('rgb(0, 0, 255)')).toBe(Color.fromRgb(0, 0, 255));
        });

        it('returns undefined for invalid colors', () => {
            expect(decodeColor('#zzzzzz')).toBeUndefined();
            expect(decodeColor('notacolor')).toBeUndefined();
            expect(decodeColor('')).toBeUndefined();
            expect(decodeColor(null)).toBeUndefined();
            expect(decodeColor(undefined)).toBeUndefined();
        });
    });

    describe('decodeColorWithAlpha', () => {
        it('6-digit hex codes with default alpha', () => {
            const result = decodeColorWithAlpha('#ff0000');
            expect(result).toBeDefined();
            expect(result?.color).toBe(Color.fromRgb(255, 0, 0));
            expect(result?.alpha).toBe(1.0);
        });

        it('8-digit hex codes with alpha', () => {
            const resultFull = decodeColorWithAlpha('#ff0000ff');
            expect(resultFull).toBeDefined();
            expect(resultFull?.color).toBe(Color.fromRgb(255, 0, 0));
            expect(resultFull?.alpha).toBeCloseTo(1.0, 2);

            const resultHalf = decodeColorWithAlpha('#00ff0080');
            expect(resultHalf).toBeDefined();
            expect(resultHalf?.color).toBe(Color.fromRgb(0, 255, 0));
            expect(resultHalf?.alpha).toBeCloseTo(0.502, 2);

            const resultZero = decodeColorWithAlpha('#0000ff00');
            expect(resultZero).toBeDefined();
            expect(resultZero?.color).toBe(Color.fromRgb(0, 0, 255));
            expect(resultZero?.alpha).toBe(0.0);

            const result7f = decodeColorWithAlpha('#ffffff7f');
            expect(result7f).toBeDefined();
            expect(result7f?.color).toBe(Color.fromRgb(255, 255, 255));
            expect(result7f?.alpha).toBeCloseTo(0.498, 2);
        });

        it('3-digit hex codes with default alpha', () => {
            const result = decodeColorWithAlpha('#f00');
            expect(result).toBeDefined();
            expect(result?.color).toBe(Color.fromRgb(255, 0, 0));
            expect(result?.alpha).toBe(1.0);
        });

        it('named colors with default alpha', () => {
            const result = decodeColorWithAlpha('red');
            expect(result).toBeDefined();
            expect(result?.color).toBe(0xff0000);
            expect(result?.alpha).toBe(1.0);
        });

        it('rgb() format with default alpha', () => {
            const result = decodeColorWithAlpha('rgb(255, 0, 0)');
            expect(result).toBeDefined();
            expect(result?.color).toBe(Color.fromRgb(255, 0, 0));
            expect(result?.alpha).toBe(1.0);
        });

        it('returns undefined for invalid colors', () => {
            expect(decodeColorWithAlpha('#zzzzzz')).toBeUndefined();
            expect(decodeColorWithAlpha('notacolor')).toBeUndefined();
            expect(decodeColorWithAlpha(null)).toBeUndefined();
            expect(decodeColorWithAlpha(undefined)).toBeUndefined();
        });

        it('handles various alpha values correctly', () => {
            const tests = [
                { hex: '#ff000000', alpha: 0.0 },
                { hex: '#ff000040', alpha: 0.251 },
                { hex: '#ff000080', alpha: 0.502 },
                { hex: '#ff0000bf', alpha: 0.749 },
                { hex: '#ff0000ff', alpha: 1.0 },
            ];

            tests.forEach(({ hex, alpha }) => {
                const result = decodeColorWithAlpha(hex);
                expect(result).toBeDefined();
                expect(result?.alpha).toBeCloseTo(alpha, 2);
            });
        });
    });
});
