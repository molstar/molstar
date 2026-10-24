/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CustomPropSymbol } from '@molstar/query-language/language/symbol';
import { Type } from '@molstar/query-language/language/type';
import { QueryRuntimeTable, QuerySymbolRuntime } from '../base.js';
import type { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import { QueryContext, Structure } from '@molstar/model/model/structure';
import { compile } from '../compiler.js';
import '../table.js'; // (registers the symbols of the queries)

function runtime(name: string) {
  return QuerySymbolRuntime.Const(CustomPropSymbol('test-table', name, Type.Num), () => 1);
}

function descriptor(name: string, ...runtimes: QuerySymbolRuntime[]): CustomPropertyDescriptor<any> {
  const symbols: any = {};
  for (const r of runtimes) symbols[r.symbol.info.name] = r;
  return { name, symbols } as CustomPropertyDescriptor<any>;
}

describe('QueryRuntimeTable counting', () => {
  it('keeps a symbol until every registration is removed', () => {
    const table = new QueryRuntimeTable();
    const r = runtime('a');

    table.addSymbol(r);
    table.addSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBe(r);

    table.removeSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBe(r);

    table.removeSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();
  });

  it('treats unknown symbol removal as a no-op', () => {
    const table = new QueryRuntimeTable();
    const r = runtime('a');

    table.removeSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();

    table.addSymbol(r);
    // a different runtime with the same id is not the registered one
    table.removeSymbol(runtime('a'));
    expect(table.getRuntime(r.symbol.id)).toBe(r);

    table.removeSymbol(r);
    table.removeSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();

    // extra removals did not leave a negative count
    table.addSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBe(r);
    table.removeSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();
  });

  it('keeps custom property symbols while another registration remains', () => {
    const table = new QueryRuntimeTable();
    const r1 = runtime('p1');
    const r2 = runtime('p2');
    const desc = descriptor('test-prop', r1, r2);

    table.addCustomProp(desc);
    table.addCustomProp(desc);
    expect(table.getRuntime(r1.symbol.id)).toBe(r1);
    expect(table.getRuntime(r2.symbol.id)).toBe(r2);

    table.removeCustomProp(desc);
    expect(table.getRuntime(r1.symbol.id)).toBe(r1);
    expect(table.getRuntime(r2.symbol.id)).toBe(r2);

    table.removeCustomProp(desc);
    expect(table.getRuntime(r1.symbol.id)).toBeUndefined();
    expect(table.getRuntime(r2.symbol.id)).toBeUndefined();
  });

  it('treats unknown custom property removal as a no-op', () => {
    const table = new QueryRuntimeTable();
    const r = runtime('p');
    const desc = descriptor('test-prop', r);

    table.removeCustomProp(desc);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();

    table.addCustomProp(desc);
    table.removeCustomProp(descriptor('test-prop', r));
    expect(table.getRuntime(r.symbol.id)).toBe(r);

    table.removeCustomProp(desc);
    table.removeCustomProp(desc);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();

    // the property can be registered again
    table.addCustomProp(desc);
    expect(table.getRuntime(r.symbol.id)).toBe(r);
  });

  it('does not remove a symbol added directly when a custom property sharing it is removed', () => {
    const table = new QueryRuntimeTable();
    const r = runtime('shared');
    const desc = descriptor('test-prop', r);

    table.addSymbol(r);
    table.addCustomProp(desc);
    table.removeCustomProp(desc);
    expect(table.getRuntime(r.symbol.id)).toBe(r);

    table.removeSymbol(r);
    expect(table.getRuntime(r.symbol.id)).toBeUndefined();
  });
});

describe('query table arguments', () => {
  it('an object-form argument list keeps its first element (perf-audit-2 06/B5)', () => {
    const run = (e: Expression) => compile(e)(new QueryContext(Structure.Empty));
    const asArray = MS.struct.type.ringFingerprint(['C', 'N', 'O']);
    const asObject: Expression = { head: (asArray as Expression.Apply).head, args: { 0: 'C', 1: 'N', 2: 'O' } };
    expect(run(asObject)).toBe(run(asArray));
    expect(run(asObject)).not.toBe(run(MS.struct.type.ringFingerprint(['N', 'O'])));
  });
});
