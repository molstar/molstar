/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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
