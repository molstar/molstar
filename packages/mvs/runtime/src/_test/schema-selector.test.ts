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
import { Column } from '@molstar/core/data/db';
import { Structure, StructureSelection, Unit } from '@molstar/model/model/structure';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';
import { SelectorParams, substructureFromSelector } from '../components/selector.js';
import { createModel } from '../helpers/_test/annotation-rows-fixture.js';
import { rowsToExpression } from '../helpers/selections.js';

setFSModule(fs);

const rowsA = [
  { label_asym_id: 'A', beg_label_seq_id: 1, end_label_seq_id: 20 },
  { label_asym_id: 'A', label_seq_id: 60, label_atom_id: 'CA' },
];
const rowsB = [{ label_comp_id: 'REA' }]; // the retinoic acid ligand of 1cbs
const rowC = { end_label_seq_id: 20 }; // the residues of all chains which have a label_seq_id up to 20 (not the ligand and the waters)
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
        {
          kind: 'component',
          params: { selector: rowC },
          children: [{ kind: 'representation', params: { type: 'ball_and_stick' } }],
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
  const components = cells.filter((c) => c.transform.transformer.id === 'ms-plugin.structure-component');
  const schemaCells = components.filter((c) => (c.transform.params as any).type.name === 'schema');
  const parentStructure = (cell: (typeof cells)[number]) =>
    plugin.state.data.cells.get(cell.transform.parent)!.obj!.data as Structure;

  it('loads without errors and creates one schema component per rows selector', () => {
    expect(cells.filter((c) => c.status === 'error').map((c) => c.errorText)).toEqual([]);
    expect(schemaCells.length).toBe(3);
    expect(schemaCells.map((c) => (c.transform.params as any).type.params)).toEqual([rowsA, rowsB, [rowC]]);
    // the labels describe the rows instead of dumping them
    expect(schemaCells.map((c) => c.obj!.label)).toEqual(['A:1-20 + A60 CA', 'REA', 'res <=20']);
    // the components of the document are all the plugin's own `StructureComponent` (so the classic UI can modify them),
    // and no rows were turned into an expression
    const types = components.map((c) => (c.transform.params as any).type.name as string);
    expect(types.length).toBe(5);
    expect(types.filter((t) => t === 'static').length).toBe(2); // (polymer, ligand)
    expect(types.filter((t) => t === 'schema').length).toBe(3);
    expect(types.includes('expression')).toBe(false);
    expect(
      cells.some((c) => c.transform.transformer.id.startsWith('mvs.') && c.transform.transformer.id.includes('schema')),
    ).toBe(false);
  });

  /** unit id -> its elements */
  const unitElements = (structure: Structure) =>
    Object.fromEntries(structure.units.map((u) => [u.id, Array.from(u.elements)]));
  const viaExpression = (cell: (typeof cells)[number], rows: any[]) =>
    StructureSelection.unionStructure(
      StructureQueryHelper.createAndRun(parentStructure(cell), rowsToExpression(rows)).selection,
    );

  it('the plugin treats them as components it can modify (the classic UI: modify by the current selection)', () => {
    expect(schemaCells.length).toBe(3);
    for (const c of schemaCells) expect(plugin.builders.structure.isComponentTransform(c)).toBe(true);
  });

  it('components contain exactly the elements the MolScript expression for their rows selects', () => {
    const [a, b] = schemaCells;
    for (const [cell, rows] of [
      [a, rowsA],
      [b, rowsB],
    ] as const) {
      const expected = viaExpression(cell, rows as any);
      expect(expected.elementCount).toBeGreaterThan(10);
      expect(unitElements(cell.obj!.data as Structure)).toEqual(unitElements(expected));
    }
  });

  it('residues without label_seq_id (the ligand, the waters) are not selected by a sequence range, unlike with MolScript', () => {
    const c = schemaCells[2];
    const structure = c.obj!.data as Structure;
    // independent of both: residues which have a label_seq_id of 20 or less, counted from the hierarchy
    const parent = parentStructure(c);
    let expected = 0;
    for (const unit of parent.units) {
      if (!Unit.isAtomic(unit)) continue;
      const { label_seq_id } = unit.model.atomicHierarchy.residues;
      for (let i = 0; i < unit.elements.length; i++) {
        const iRes = unit.residueIndex[unit.elements[i]];
        if (label_seq_id.valueKind(iRes) === Column.ValueKind.Present && label_seq_id.value(iRes) <= 20) expected++;
      }
    }
    expect(expected).toBeGreaterThan(100);
    expect(structure.elementCount).toBe(expected);
    expect(viaExpression(c, [rowC]).elementCount).toBeGreaterThan(expected); // (MolScript reads a missing label_seq_id as 0)
  });

  it('every schema component has a representation, which was built', () => {
    const descendantsOf = (ref: string): typeof cells => {
      const result: typeof cells = [];
      const queue = [ref];
      while (queue.length) {
        const parent = queue.pop()!;
        for (const c of cells) {
          if (c.transform.parent === parent) {
            result.push(c);
            queue.push(c.transform.ref);
          }
        }
      }
      return result;
    };
    for (const c of schemaCells) {
      const representations = descendantsOf(c.transform.ref).filter(
        (d) => d.transform.transformer.id === 'ms-plugin.structure-representation-3d',
      );
      expect(representations.length).toBeGreaterThan(0);
      expect(representations.map((r) => r.status)).toEqual(representations.map(() => 'ok'));
    }
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
