/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import * as fs from 'node:fs';
import { setFSModule } from '@molstar/core/util/data-source';
import { StateObject, StateTransformer } from '@molstar/core/state';
import { PluginContext } from '@molstar/plugin/context';
import { PluginSpec } from '@molstar/plugin/spec';
import { Mmcif } from '@molstar/plugin/state/formats/trajectory/mmcif';
import { StructureQueryHelper } from '@molstar/plugin/state/helpers/structure-query';
import { Structure, StructureSelection } from '@molstar/model/model/structure';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';
import {
  createMVSSchemaStructureComponent,
  updateMVSSchemaStructureComponent,
} from '../components/schema-structure-component.js';
import { SelectorParams, substructureFromSelector } from '../components/selector.js';
import { createModel } from '../helpers/_test/annotation-rows-fixture.js';
import { rowsToExpression } from '../helpers/selections.js';

setFSModule(fs);

const rowsA = [
  { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 20 },
  { label_asym_id: 'A', label_seq_id: 60, label_atom_id: 'CA' },
];
const rowsB = [{ label_comp_id: 'REA' }]; // the retinoic acid ligand of 1cbs
const colorLayerRows = { label_asym_id: 'A', beg_label_seq_id: 21, end_label_seq_id: 40 };

/** `data/examples/mvs/1cbs.mvsj` (local structure file, no labels: Node has no canvas) with extra components selected by rows. */
function createDocument() {
  const mvsj = JSON.parse(fs.readFileSync('data/examples/mvs/1cbs.mvsj', 'utf8'));
  const visit = (node: any) => {
    if (node.kind === 'download') node.params.url = `file://${process.cwd()}/data/examples/1cbs_full.bcif`;
    node.children = node.children?.filter((c: any) => c.kind !== 'label');
    node.children?.forEach(visit);
    if (node.kind === 'structure') {
      node.children.push(
        {
          kind: 'component',
          params: { selector: rowsA },
          children: [
            { kind: 'representation', params: { type: 'ball_and_stick' } },
            { kind: 'tooltip', params: { text: 'rows A' } },
          ],
        },
        {
          kind: 'component',
          params: { selector: rowsB[0] },
          children: [{ kind: 'representation', params: { type: 'spacefill' } }],
        },
      );
    }
  };
  visit(mvsj.root);
  // a color layer with a row selector on the existing cartoon representation
  const repr = mvsj.root.children[0].children[0].children[0].children[0].children[0];
  repr.children.push({ kind: 'color', params: { selector: colorLayerRows, color: 'red' } });
  return MVSData.fromMVSJ(JSON.stringify(mvsj));
}

describe('rows selectors in a loaded MVS document', async () => {
  const plugin = new PluginContext({
    registry: [MVSRuntimeRegistry, Mmcif],
    behaviors: [PluginSpec.Behavior(MolViewSpec)],
  });
  await plugin.init();
  const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
  const error = vi.spyOn(console, 'error').mockImplementation(() => {});
  try {
    await loadMVS(plugin, createDocument(), { sanityChecks: true, keepCamera: true });
  } finally {
    warn.mockRestore();
    error.mockRestore();
  }
  const cells = [...plugin.state.data.cells.values()];
  const schemaCells = cells.filter((c) => c.transform.transformer.id === 'mvs.mvs-structure-component-from-schema');
  const parentStructure = (cell: (typeof cells)[number]) =>
    plugin.state.data.cells.get(cell.transform.parent)!.obj!.data as Structure;
  const referenceElements = (cell: (typeof cells)[number], rows: any[]) =>
    StructureSelection.unionStructure(
      StructureQueryHelper.createAndRun(parentStructure(cell), rowsToExpression(rows)).selection,
    ).elementCount;

  it('loads without errors and creates one schema component per rows selector', () => {
    expect(cells.filter((c) => c.status === 'error').map((c) => c.errorText)).toEqual([]);
    expect(schemaCells.length).toBe(2);
    expect(schemaCells.map((c) => (c.transform.params as any).schema)).toEqual([rowsA, rowsB]);
    // plugin's own `StructureComponent` is only used for static selectors, rows are not turned into expressions
    const standard = cells.filter((c) => c.transform.transformer.id === 'ms-plugin.structure-component');
    expect(standard.map((c) => (c.transform.params as any).type.name)).toEqual(['static', 'static']);
  });

  it('components contain exactly the elements the MolScript expression for their rows selects', () => {
    const [a, b] = schemaCells;
    const expectedA = referenceElements(a, rowsA);
    const expectedB = referenceElements(b, rowsB);
    expect(expectedA).toBeGreaterThan(100);
    expect(expectedB).toBeGreaterThan(10);
    expect((a.obj!.data as Structure).elementCount).toBe(expectedA);
    expect((b.obj!.data as Structure).elementCount).toBe(expectedB);
  });

  it('the representations on schema components were built', () => {
    for (const c of schemaCells) {
      const reprs = plugin.state.data.select(c.transform.ref).length; // the component itself resolves
      expect(reprs).toBe(1);
    }
    const reprCells = cells.filter((c) => c.transform.transformer.id === 'ms-plugin.structure-representation-3d');
    expect(reprCells.every((c) => c.status === 'ok')).toBe(true);
  });

  it('tooltips and color layers carry the rows as `schema` selectors', () => {
    const props = cells.find((c) => c.transform.transformer.id === 'ms-plugin.custom-structure-properties')!;
    const tooltips = JSON.stringify((props.transform.params as any).properties);
    expect(tooltips).toContain('"name":"schema"');
    expect(tooltips).toContain('rows A');
    const cartoon = cells.find((c) => (c.transform.params as any)?.type?.name === 'cartoon')!;
    const layers = (cartoon.transform.params as any).colorTheme.params.layers as { selection: any }[];
    expect(layers.map((l) => l.selection.name)).toContain('schema');
    expect(layers.map((l) => l.selection.params)).toContainEqual([colorLayerRows]);
  });
});

describe('MVSSchemaStructureComponent', async () => {
  const dims = { nChains: 4, nResidues: 20 };
  const structure = Structure.ofModel(await createModel(dims));
  const rows = [{ label_asym_id: 'B' }];

  it('creates a component with the rows` elements, a label, and `null` when empty and nullIfEmpty', () => {
    const cache: { source?: Structure } = {};
    const obj = createMVSSchemaStructureComponent(
      structure,
      { schema: rows, label: '', nullIfEmpty: true },
      cache,
    ) as any;
    expect(obj.data.elementCount).toBe(20 * 6);
    expect(obj.label).toBe('MVS Schema Component');
    expect(cache.source).toBe(structure);
    const named = createMVSSchemaStructureComponent(
      structure,
      { schema: rows, label: 'Chain B', nullIfEmpty: true },
      {},
    ) as any;
    expect(named.label).toBe('Chain B');
    const none = { schema: [{ label_asym_id: 'nonexistent' }], label: '', nullIfEmpty: true };
    expect(createMVSSchemaStructureComponent(structure, none, {})).toBe(StateObject.Null);
    expect(
      (createMVSSchemaStructureComponent(structure, { ...none, nullIfEmpty: false }, {}) as any).data.elementCount,
    ).toBe(0);
  });

  it('update: label only -> Updated, same params -> Unchanged, other rows or other structure -> Recreate', () => {
    const params = { schema: rows, label: 'x', nullIfEmpty: true };
    const cache: { source?: Structure } = {};
    const obj = createMVSSchemaStructureComponent(structure, params, cache) as any;
    const U = StateTransformer.UpdateResult;
    expect(updateMVSSchemaStructureComponent(structure, obj, params, { ...params }, cache)).toBe(U.Unchanged);
    expect(updateMVSSchemaStructureComponent(structure, obj, params, { ...params, label: 'y' }, cache)).toBe(U.Updated);
    expect(obj.label).toBe('y');
    expect(
      updateMVSSchemaStructureComponent(structure, obj, params, { ...params, schema: [{ label_asym_id: 'C' }] }, cache),
    ).toBe(U.Recreate);
    const smaller = substructureFromSelector(structure, { name: 'schema', params: rows });
    expect(updateMVSSchemaStructureComponent(smaller, obj, params, params, cache)).toBe(U.Recreate);
  });
});

describe('Selector with schema rows', async () => {
  const structure = Structure.ofModel(await createModel({ nChains: 4, nResidues: 20 }));

  it('is a valid `SelectorParams` value, with an empty default', () => {
    expect(SelectorParams.select.options.map((o) => o[0])).toContain('schema');
    expect(SelectorParams.map('schema').defaultValue).toEqual([]);
  });
  it('substructureFromSelector resolves rows like the expression selector does', () => {
    const rows = [{ label_asym_id: 'A', beg_label_seq_id: 3, end_label_seq_id: 9 }, { auth_asym_id: 'x2' }];
    const viaSchema = substructureFromSelector(structure, { name: 'schema', params: rows });
    const viaExpression = substructureFromSelector(structure, { name: 'expression', params: rowsToExpression(rows) });
    expect(viaSchema.elementCount).toBe(viaExpression.elementCount);
    expect(viaSchema.elementCount).toBeGreaterThan(0);
  });
});
