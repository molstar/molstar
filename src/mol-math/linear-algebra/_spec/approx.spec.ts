/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { fasterPow2, fasterExp, fasterLog, fasterLog10, fasterSin, fasterCos, fastAtan, fastAtan2, fasterTan, fasterTanh, fasterCosh, fasterSinh, fastPow2, fastExp, fastLog, fastLog10, fastSinh, fastCosh, fastTanh, fastSin, fastCos, fastTan } from '../../approx';

describe('approx', () => {
    it('fastPow2', () => {
        expect(fastPow2(4)).toBeCloseTo(Math.pow(2, 4), 2);
    });

    it('fasterPow2', () => {
        expect(fasterPow2(4)).toBeCloseTo(Math.pow(2, 4), 0);
    });

    it('fastExp', () => {
        expect(fastExp(4)).toBeCloseTo(Math.exp(4), 2);
    });

    it('fasterExp', () => {
        expect(fasterExp(4)).toBeCloseTo(Math.exp(4), 0);
    });

    it('fastLog', () => {
        expect(fastLog(12)).toBeCloseTo(Math.log(12), 2);
    });

    it('fasterLog', () => {
        expect(fasterLog(12)).toBeCloseTo(Math.log(12), 1);
    });

    it('fastLog10', () => {
        expect(fastLog10(42)).toBeCloseTo(Math.log10(42), 2);
    });

    it('fasterLog10', () => {
        expect(fasterLog10(42)).toBeCloseTo(Math.log10(42), 1);
    });

    it('fastSinh', () => {
        expect(fastSinh(0.3)).toBeCloseTo(Math.sinh(0.3), 2);
    });

    it('fasterSinh', () => {
        expect(fasterSinh(0.3)).toBeCloseTo(Math.sinh(0.3), 1);
    });

    it('fastCosh', () => {
        expect(fastCosh(0.3)).toBeCloseTo(Math.cosh(0.3), 2);
    });

    it('fasterCosh', () => {
        expect(fasterCosh(0.3)).toBeCloseTo(Math.cosh(0.3), 1);
    });

    it('fastTanh', () => {
        expect(fastTanh(0.3)).toBeCloseTo(Math.tanh(0.3), 2);
    });

    it('fasterTanh', () => {
        expect(fasterTanh(0.3)).toBeCloseTo(Math.tanh(0.3), 1);
    });

    it('fastSin', () => {
        expect(fastSin(0.3)).toBeCloseTo(Math.sin(0.3), 2);
    });

    it('fasterSin', () => {
        expect(fasterSin(0.3)).toBeCloseTo(Math.sin(0.3), 1);
    });

    it('fastCos', () => {
        expect(fastCos(0.3)).toBeCloseTo(Math.cos(0.3), 2);
    });

    it('fasterCos', () => {
        expect(fasterCos(0.3)).toBeCloseTo(Math.cos(0.3), 1);
    });

    it('fastTan', () => {
        expect(fastTan(0.3)).toBeCloseTo(Math.tan(0.3), 2);
    });

    it('fasterTan', () => {
        expect(fasterTan(0.3)).toBeCloseTo(Math.tan(0.3), 1);
    });

    it('fastAtan', () => {
        expect(fastAtan(0.3)).toBeCloseTo(Math.atan(0.3), 2);
    });

    it('fastAtan2', () => {
        expect(fastAtan2(0.1, 0.5)).toBeCloseTo(Math.atan2(0.1, 0.5), 2);
    });
});