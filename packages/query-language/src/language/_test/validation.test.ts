/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { MolScriptBuilder as B } from '../builder.js';
import { Expression } from '../expression.js';
import { CustomPropSymbol } from '../symbol.js';
import { SymbolMap } from '../symbol-table.js';
import { Type } from '../type.js';
import { expressionValidationIssues as issues } from '../validation.js';

describe('MolQL syntax and field names', () => {
  it('accepts nested queries, positional arrays/maps, and bare-symbol string values', () => {
    const expression = B.struct.generator.atomGroups({
      'chain-test': B.core.rel.eq({ 0: B.ammp('label_asym_id'), 1: Expression.Symbol('A') }),
      'atom-test': B.core.logic.and([B.core.rel.eq([B.acp('elementSymbol'), B.es('C')])]),
    });
    expect(issues(expression)).toBeUndefined();
    expect(issues(B.re('^C'))).toBeUndefined();
    expect(issues(JSON.parse(JSON.stringify(B.re('^C'))))).toBeUndefined();
  });

  it.each([
    null,
    [],
    {},
    { head: null },
    { head: 'call' },
    { head: { name: 1 } },
    { head: { name: 'core.math.add' }, args: 1 },
  ])('rejects malformed syntax %p', (expression) => expect(issues(expression)).toBeDefined());

  it('rejects unexpected expression fields even in syntax-only mode', () => {
    expect(
      issues({ head: { name: 'custom.call', extra: true }, args: [], arguments: [] }, { syntaxOnly: true })?.join('\n'),
    ).toContain("Unknown application field 'arguments'");
    expect(issues({ name: 'A', args: [] })?.join('\n')).toContain("Unknown symbol field 'args'");
    expect(issues({ head: { name: 'core.math.add' }, args: [Infinity] })?.join('\n')).toContain('Expected a literal');
  });

  it('accepts null and lists of values as arguments (the columns of a table), not as an expression of their own', () => {
    const table = B.struct.type.schemaItem({
      label_asym_id: ['A', null, 'B'],
      label_seq_id: [1, 2, null],
      auth_seq_id: null,
    });
    expect(issues(B.struct.generator.elementSchema([table]))).toBeUndefined();
    expect(issues(JSON.parse(JSON.stringify(B.struct.generator.elementSchema([table]))))).toBeUndefined();
    expect(issues(B.struct.type.schemaItem({ label_seq_id: [1, Infinity] }))?.join('\n')).toContain(
      'expression.args["label_seq_id"][1]: Expected a string, a finite number, a boolean or null in a list of values.',
    );
    expect(issues(B.struct.type.schemaItem({ label_seq_id: [1, [2]] as any }))).toBeDefined();
    expect(issues(B.struct.type.schemaItem({ label_seq_id: [1, { name: 'A' }] as any }))).toBeDefined();
    expect(issues([1, 2])).toBeDefined();
    expect(issues(null)).toBeDefined();
  });

  it('accepts a list of values or null only where the argument can be one', () => {
    // (a list or a missing value is not a boolean)
    expect(issues(B.struct.generator.atomGroups({ 'residue-test': [false] as any }))?.join('\n')).toContain(
      'expression.args["residue-test"]: This argument is not a list of values, or a missing value (null).',
    );
    expect(issues(B.struct.generator.atomGroups({ 'chain-test': null as any }))).toBeDefined();
    expect(issues(B.core.math.add([[1, 2], 1] as any))).toBeDefined();
    // (a list is what a list, a set and a field of an item which can be a column are)
    expect(issues(B.core.type.list([[1, 2], null] as any))).toBeUndefined();
    expect(issues(B.core.type.set([[1, 2]] as any))).toBeUndefined();
    expect(issues(B.struct.type.schemaItem({ label_seq_id: [1, null], label_asym_id: null }))).toBeUndefined();
    // (and of any type)
    expect(issues(B.struct.generator.atomGroups({ 'group-by': [1, 2] as any }))).toBeUndefined();
    // (nothing is known about a symbol which is not looked up)
    expect(issues({ head: { name: 'custom.call' }, args: { x: [1, null] } }, { syntaxOnly: true })).toBeUndefined();
    // (and the elements of the generator are items, not lists)
    expect(issues(B.struct.generator.elementSchema([[1, 2]] as any))).toBeDefined();
  });

  it('rejects unknown nested callable symbols with their paths', () => {
    const expression = B.struct.generator.atomGroups({ 'atom-test': { head: { name: 'unknown' } } });
    expect(issues(expression)?.join('\n')).toContain('expression.args["atom-test"].head.name');
    expect(issues(expression)?.join('\n')).toContain("Unknown callable symbol 'unknown'");
    expect(issues(expression, { syntaxOnly: true })).toBeUndefined();
  });

  it('rejects misspelled and undeclared argument keys', () => {
    expect(issues(B.struct.generator.atomGroups({ 'atom-tset': true } as any))?.join('\n')).toContain(
      "Unknown argument 'atom-tset'",
    );
    expect(issues(B.core.rel.eq([1, 2, 3]))?.join('\n')).toContain("Unknown argument '2'");
    expect(issues(B.core.math.add({ typo: 1 }))?.join('\n')).toContain("Unknown argument 'typo'");
    expect(issues(B.core.math.add({ 0: 1, 1: 2 }))).toBeUndefined();
  });

  it('does not check required arguments, value types, result types, or runtime support', () => {
    expect(issues(B.struct.filter.pick({}))).toBeUndefined();
    expect(issues(B.struct.generator.atomGroups({ 'atom-test': 'a string' }))).toBeUndefined();
    expect(issues(B.core.math.add([1, 2]))).toBeUndefined();
    expect(issues(B.core.ctrl.fn([1]))).toBeUndefined();
  });

  it('accepts explicitly supplied custom symbol definitions', () => {
    const custom = CustomPropSymbol('custom', 'value', Type.Num);
    const expression = B.core.rel.eq([custom(), 1]);
    expect(issues(expression)).toBeDefined();
    expect(
      issues(expression, { getSymbol: (name) => (name === custom.id ? custom : SymbolMap[name]) }),
    ).toBeUndefined();
  });

  it('rejects cycles without rejecting shared expression objects', () => {
    const expression: any = B.core.math.add([]);
    expression.args.push(expression);
    expect(issues(expression)?.join('\n')).toContain('Cyclic expression');
    const shared = B.ammp('label_asym_id');
    expect(issues(B.core.rel.eq([shared, shared]))).toBeUndefined();
  });

  it('rejects sparse argument arrays and properties lost during JSON serialization', () => {
    expect(issues(B.core.math.add(new Array(1)))?.join('\n')).toContain('expression.args["0"]');
    const args: any = [1];
    args.typo = 2;
    expect(issues(B.core.math.add(args), { syntaxOnly: true })?.join('\n')).toContain("Unknown array field 'typo'");
  });
});
