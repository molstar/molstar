/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Expression } from '../expression.js';

describe('MolScript expression shape', () => {
  it('accepts literals, symbols, and recursive applications', () => {
    expect(Expression.is('text')).toEqual(true);
    expect(Expression.is(42)).toEqual(true);
    expect(Expression.is(false)).toEqual(true);
    expect(Expression.is({ name: 'structure-query.generator.all' })).toEqual(true);
    expect(
      Expression.is({
        head: { name: 'core.rel.eq' },
        args: [{ name: 'structure.atom-property.macromolecular.label_asym_id' }, 'A'],
      }),
    ).toEqual(true);
    expect(
      Expression.is({
        head: { name: 'structure-query.generator.atom-groups' },
        args: { 'chain-test': { head: { name: 'core.rel.eq' }, args: ['A', 'A'] } },
      }),
    ).toEqual(true);
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

  it('null and lists of values are literals in the arguments of an application', () => {
    expect(Expression.is({ head: { name: 'a' }, args: { x: [1, null, 'b', true], y: null } })).toEqual(true);
    expect(Expression.is({ head: { name: 'a' }, args: [[1, 2], null] })).toEqual(true);
    expect(Expression.is({ head: { name: 'a' }, args: { x: [1, [2]] } })).toEqual(false);
    expect(Expression.is({ head: { name: 'a' }, args: { x: [{ name: 's' }] } })).toEqual(false);
    expect(Expression.isLiteral([1, null] as any)).toEqual(true);
    expect(Expression.isApply([1, null] as any)).toEqual(false);
    expect(Expression.isSymbol([1, null] as any)).toEqual(false);
    // (on their own they are not expressions)
    expect(Expression.is([1, 2])).toEqual(false);
    expect(Expression.is(null)).toEqual(false);
  });
});
