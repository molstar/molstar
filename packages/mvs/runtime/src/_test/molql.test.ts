/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { componentPropsFromSelector, prettyNameFromSelector } from '@molstar/mvs/load-helpers';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { createMVSBuilder } from '@molstar/mvs-builder/tree/mvs/mvs-builder';
import {
  MolQLExpressionT,
  PrimitiveMolQLExpressionT,
  isMolQLExpression,
  isPrimitiveMolQLExpression,
} from '@molstar/mvs-builder/tree/mvs/param-types';
import { molQLValidationIssues } from '@molstar/mvs/runtime-validation';
import { compile, DefaultQueryRuntimeTable, QuerySymbolRuntime } from '@molstar/model/script/runtime/query/base';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import { CustomPropSymbol } from '@molstar/query-language/language/symbol';
import { Type } from '@molstar/query-language/language/type';
import { parse } from '@molstar/query-language/transpile';
import '@molstar/query-language/transpilers/pymol';

describe('MVS MolQL selectors', () => {
  const molql = MS.struct.generator.atomGroups({});

  it('ignores nested extra fields during both name validation and executable compilation', () => {
    // The first extra passes language validation but fails executable compilation; the second has an unknown name.
    const extras = [MS.core.ctrl.fn([1]), { head: { name: 'custom.future-query' } }];
    for (const expression of extras) {
      const builder = createMVSBuilder();
      const structure = builder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure();
      const selector = { label_asym_id: 'A', extension_data: { molql: expression } };
      structure.component({ selector });
      structure.component({ selector: [selector] });
      structure.primitives().label({
        position: { expressions: [selector], extension_data: { molql: expression } } as any,
        text: 'label',
      });
      expect(molQLValidationIssues(builder.getState().root)).toBeUndefined();
    }
  });

  it('still compiles actual selectors and primitive positions', () => {
    const builder = createMVSBuilder();
    const structure = builder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure();
    const unsupported = MS.core.ctrl.fn([1]);
    structure.component({ selector: { molql: unsupported } });
    structure.primitives().label({ position: { molql: unsupported, structure_ref: 's' }, text: 'label' });
    const issues = molQLValidationIssues(builder.getState().root)?.join('\n');
    expect(issues).toContain('.params.selector.molql: Invalid MolQL expression');
    expect(issues).toContain('.params.position.molql: Invalid MolQL expression');
  });

  it('validates the standalone MolQL wrapper structure', () => {
    const malformed = { head: { name: 'not-implemented' } };
    expect(MolQLExpressionT.decode({ molql })._tag).toEqual('Right');
    expect(MolQLExpressionT.decode({ molql: 'not-a-query' })._tag).toEqual('Left');
    // Fixed schema codecs check syntax; MVSData validates names using the caller's symbol vocabulary.
    expect(MolQLExpressionT.decode({ molql: malformed })._tag).toEqual('Right');
    expect(MolQLExpressionT.decode({ molql, label_seq_id: 1 })._tag).toEqual('Left');
    expect(MolQLExpressionT.decode({})._tag).toEqual('Left');
    expect(isMolQLExpression({ molql })).toEqual(true);
    expect(isMolQLExpression({})).toEqual(false);
    expect(isMolQLExpression([])).toEqual(false);
  });

  it('allows primitive MolQL positions to reference another structure', () => {
    const position = { molql, structure_ref: 'other-structure' };
    expect(PrimitiveMolQLExpressionT.decode(position)._tag).toEqual('Right');
    expect(PrimitiveMolQLExpressionT.decode({ molql, structure_ref: 1 })._tag).toEqual('Left');
    expect(isPrimitiveMolQLExpression(position)).toEqual(true);
    expect(isMolQLExpression({ molql, structure_ref: 1 })).toEqual(false);
    expect(isPrimitiveMolQLExpression({ molql, structure_ref: 1 })).toEqual(false);
  });

  it('passes wrapped expressions through unchanged and preserves legacy selectors', () => {
    const wrapped = { molql };
    expect(componentPropsFromSelector(wrapped)).toEqual({ name: 'expression', params: molql });
    expect(componentPropsFromSelector('protein')).toEqual({ name: 'static', params: 'protein' });
    // rows are kept as rows (resolved through model indices), not converted to a MolScript expression
    expect(componentPropsFromSelector({ label_asym_id: 'A' })).toEqual({
      name: 'schema',
      params: [{ label_asym_id: 'A' }],
    });
    const rows = [{ label_asym_id: 'A' }, { label_asym_id: 'B', label_seq_id: 3 }];
    expect(componentPropsFromSelector(rows)).toEqual({ name: 'schema', params: rows });
    expect(prettyNameFromSelector(wrapped)).toEqual('MolQL Selection');
    expect(prettyNameFromSelector({ label_asym_id: 'A' })).toEqual('Chain A');
    expect(prettyNameFromSelector(rows)).toEqual('Chain A + Chain B, residue 3');
    expect(prettyNameFromSelector('protein')).toEqual('Protein');
    expect(prettyNameFromSelector(undefined)).toEqual('All');
  });

  it('rejects invalid selector and primitive addresses during MVS validation', () => {
    const builder = createMVSBuilder();
    builder
      .download({ url: 'example.bcif' })
      .parse({ format: 'bcif' })
      .modelStructure()
      .component({ selector: { molql: 'not-a-query' } as any });
    expect(MVSData.validationIssues(builder.getState())?.join('\n')).toContain(
      'MolQL expression must be a MolScript application',
    );

    const mixedBuilder = createMVSBuilder();
    mixedBuilder
      .download({ url: 'example.bcif' })
      .parse({ format: 'bcif' })
      .modelStructure()
      .component({ selector: { molql, label_seq_id: 1 } as any });
    expect(MVSData.validationIssues(mixedBuilder.getState())?.join('\n')).toContain(
      'Expected an object with only a molql property',
    );

    const primitiveBuilder = createMVSBuilder();
    const structure = primitiveBuilder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure();
    structure.primitives().label({ position: { molql, structure_ref: 1 } as any, text: 'Invalid reference' });
    expect(MVSData.validationIssues(primitiveBuilder.getState())?.join('\n')).toContain(
      'Expected an object with only molql and an optional string structure_ref',
    );

    const referencedBuilder = createMVSBuilder();
    const referenced = referencedBuilder
      .download({ url: 'other.bcif' })
      .parse({ format: 'bcif' })
      .modelStructure({ ref: 'struct1' });
    referenced
      .primitives()
      .label({ position: { structure_ref: 'struct1', label_asym_id: 'A' }, text: 'Legacy reference' });
    expect(MVSData.validationIssues(referencedBuilder.getState())).toEqual(undefined);
  });

  it('retains MolScript compiler validation in the runtime layer', () => {
    const builder = createMVSBuilder();
    builder
      .download({ url: 'example.bcif' })
      .parse({ format: 'bcif' })
      .modelStructure()
      .component({ selector: { molql: { head: { name: 'not-implemented' } } } as any });
    const state = builder.getState();
    expect(MVSData.validationIssues(state)?.join('\n')).toContain('Unknown callable symbol');
    expect(molQLValidationIssues(state.root)?.join('\n')).toContain('Invalid MolQL expression');
  });

  it('validates registered custom symbols in the runtime without changing standalone vocabulary', () => {
    const symbol = CustomPropSymbol('custom', 'mvs-test', Type.Num);
    const runtime = QuerySymbolRuntime.Const(symbol, () => 1);
    const builder = createMVSBuilder();
    builder
      .download({ url: 'example.bcif' })
      .parse({ format: 'bcif' })
      .modelStructure()
      .component({
        selector: { molql: MS.struct.generator.atomGroups({ 'atom-test': MS.core.rel.eq([symbol(), 1]) }) },
      });
    const state = builder.getState();
    DefaultQueryRuntimeTable.addSymbol(runtime);
    try {
      expect(molQLValidationIssues(state.root)).toBeUndefined();
      expect(MVSData.validationIssues(state)?.join('\n')).toContain('Unknown callable symbol');
    } finally {
      DefaultQueryRuntimeTable.removeSymbol(runtime);
    }
    expect(molQLValidationIssues(state.root)?.join('\n')).toContain('Unknown callable symbol');
  });

  it('compiles MolScriptBuilder and PyMOL-transpiled expressions', () => {
    expect(() => compile(parse('pymol', 'byres polymer within 5 of resn STI'))).not.toThrow();
    expect(() =>
      compile(
        MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), 'G']),
          'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), 'N13']),
        }),
      ),
    ).not.toThrow();
    expect(() => compile(parse('pymol', 'chain A and resi 315 and name OG1'))).not.toThrow();
  });
});
