/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData } from '../mvs-data.js';
import { MolScriptBuilder as B, compileScript } from '../molql.js';
import { molQLValidationIssues } from '../molql-validation.js';
import { MolQLExpressionT } from '../tree/mvs/param-types.js';
import { CustomPropSymbol } from '@molstar/query-language/language/symbol';
import { SymbolMap } from '@molstar/query-language/language/symbol-table';
import { Type } from '@molstar/query-language/language/type';
import type { Expression } from '@molstar/query-language/language/expression';

function state(expression: Expression) {
  const builder = MVSData.createBuilder();
  builder
    .download({ url: 'example.bcif' })
    .parse({ format: 'bcif' })
    .modelStructure()
    .component({ selector: { molql: expression } });
  return builder.getState();
}

describe('Standalone MVS MolQL authoring and validation', () => {
  it.each(['mol-script', 'pymol', 'vmd', 'jmol'] as const)('builds and serializes %s selectors', (language) => {
    const expression = compileScript(language, language === 'mol-script' ? '(sel.atom.all)' : 'all');
    const data = MVSData.fromMVSJ(MVSData.toMVSJ(state(expression)));
    expect(MVSData.isValid(data)).toBe(true);
    expect(MVSData.validationIssues(data)).toBeUndefined();
  });

  it('rejects unknown calls and misspelled argument names in standalone validation', () => {
    expect(MVSData.validationIssues(state({ head: { name: 'unknown' } }))?.join('\n')).toContain(
      'Unknown callable symbol',
    );
    const data = state(B.struct.generator.atomGroups({ 'chain-tset': true } as any));
    expect(MVSData.isValid(data)).toBe(false);
    expect(MVSData.validationIssues(data)?.join('\n')).toContain("Unknown argument 'chain-tset'");
  });

  it('ignores MolQL-looking values in extra parameters unless extra parameters are forbidden', () => {
    const data = state(B.struct.generator.all());
    const extra = { molql: { head: { name: 'custom.future-query' } } };
    (data.root as any).params = { extension_data: extra };
    expect(MVSData.validationIssues(data)).toBeUndefined();
    expect(MVSData.validationIssues(data, { noExtra: false })).toBeUndefined();
    expect(MVSData.validationIssues(data, { noExtra: true })?.join('\n')).toContain(
      'Unknown parameter "extension_data"',
    );

    const component = data.root.children![0].children![0].children![0].children![0];
    (component.params as any).selector = extra;
    expect(MVSData.validationIssues(data)?.join('\n')).toContain("Unknown callable symbol 'custom.future-query'");
  });

  it('ignores extra parameters in discriminated primitive schemas and snapshots', () => {
    const builder = MVSData.createBuilder();
    const structure = builder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure();
    structure.primitives().label({ position: [0, 0, 0], text: 'label' });
    const data = MVSData.stateToStates(builder.getState());
    const primitive = data.snapshots[0].root.children![0].children![0].children![0].children![0].children![0];
    (primitive.params as any).extension_data = { molql: { head: { name: 'unknown' } } };
    expect(MVSData.validationIssues(data)).toBeUndefined();
    expect(MVSData.validationIssues(data, { noExtra: true })?.join('\n')).toContain(
      'Unknown parameter "extension_data"',
    );
    (primitive.params as any).position = { molql: { head: { name: 'unknown' } } };
    expect(MVSData.validationIssues(data)?.join('\n')).toContain('.position.molql.head.name');
  });

  it('validates every snapshot and primitive position with its path', () => {
    const builder = MVSData.createBuilder();
    const structure = builder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure({ ref: 's' });
    structure
      .primitives()
      .label({ position: { molql: { head: { name: 'unknown' } }, structure_ref: 's' }, text: 'label' });
    const data = MVSData.stateToStates(builder.getState());
    expect(MVSData.validationIssues(data)?.join('\n')).toContain('snapshots[0].root.children[0]');
    expect(MVSData.validationIssues(data)?.join('\n')).toContain('.position.molql.head.name');
  });

  it('checks expression wrappers nested in animation parameters', () => {
    const tree = { kind: 'animation', params: { end: [{ molql: { head: { name: 'unknown' } } }] } };
    expect(molQLValidationIssues(tree, {}, 'animation')?.join('\n')).toContain(
      'animation.params.end[0].molql.head.name',
    );
  });

  it('allows a custom vocabulary without changing fixed schema codecs', () => {
    const symbol = CustomPropSymbol('custom', 'value', Type.Num);
    const expression = B.struct.generator.atomGroups({ 'atom-test': B.core.rel.eq([symbol(), 1]) });
    expect(MolQLExpressionT.decode({ molql: expression })._tag).toBe('Right');
    const data = state(expression);
    expect(MVSData.isValid(data)).toBe(false);
    expect(MVSData.isValid(data, { getMolQLSymbol: (name) => (name === symbol.id ? symbol : SymbolMap[name]) })).toBe(
      true,
    );
  });

  it('reports malformed snapshot containers as issues', () => {
    expect(
      MVSData.validationIssues({ kind: 'multiple', metadata: { version: '1' }, snapshots: null } as any),
    ).toBeDefined();
    expect(
      MVSData.validationIssues({ kind: 'multiple', metadata: { version: '1' }, snapshots: [null] } as any),
    ).toBeDefined();
  });

  it('rejects cyclic expressions through the schema codec without overflowing', () => {
    const expression: any = B.core.math.add([]);
    expression.args.push(expression);
    expect(MolQLExpressionT.decode({ molql: expression })._tag).toBe('Left');
  });
});
