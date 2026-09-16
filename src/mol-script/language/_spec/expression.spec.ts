/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Expression } from '../expression';

describe('MolScript expression shape', () => {
    it('accepts literals, symbols, and recursive applications', () => {
        expect(Expression.is('text')).toEqual(true);
        expect(Expression.is(42)).toEqual(true);
        expect(Expression.is(false)).toEqual(true);
        expect(Expression.is({ name: 'structure-query.generator.all' })).toEqual(true);
        expect(Expression.is({
            head: { name: 'core.rel.eq' },
            args: [{ name: 'structure.atom-property.macromolecular.label_asym_id' }, 'A'],
        })).toEqual(true);
        expect(Expression.is({
            head: { name: 'structure-query.generator.atom-groups' },
            args: { 'chain-test': { head: { name: 'core.rel.eq' }, args: ['A', 'A'] } },
        })).toEqual(true);
    });

    it('rejects malformed recursive shapes', () => {
        expect(Expression.is(null)).toEqual(false);
        expect(Expression.is([])).toEqual(false);
        expect(Expression.is({})).toEqual(false);
        expect(Expression.is({ head: null })).toEqual(false);
        expect(Expression.is({ head: { name: 'core.rel.eq' }, args: [undefined] })).toEqual(false);
        expect(Expression.is({ head: { name: 'core.rel.eq' }, args: 1 })).toEqual(false);
        expect(Expression.is({ head: { malformed: true } })).toEqual(false);
    });
});
