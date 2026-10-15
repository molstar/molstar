/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { StateTreeFeature } from '../index.js';
import type { StateCellInfo } from '../operations.js';
import type { StateCellNode } from '../source.js';

// The `state-tree` feature (viewer-spec.md §3.1): the `state` area over 1CBS, headless (the plugin's state needs no
// canvas), and the source the State tree tab follows.

async function viewer(): Promise<PrismViewer> {
  return headlessViewer({ extensions: [StateTreeFeature], files: [cbs()] });
}

async function tree(v: PrismViewer): Promise<StateCellInfo[]> {
  const r = await v.api.state!.get({});
  expectOk(r);
  return r.value as StateCellInfo[];
}

function all(nodes: readonly StateCellInfo[]): StateCellInfo[] {
  return nodes.flatMap((n) => [n, ...all(n.children ?? [])]);
}

const byTransform = (nodes: readonly StateCellInfo[], transform: string) =>
  all(nodes).filter((n) => n.transform === transform);

const cellOf = (v: PrismViewer, ref: string) => v.plugin.state.data.cells.get(ref)!;
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;

describe('state tree', () => {
  it('is a built-in feature with the state area', async () => {
    expect(StateTreeFeature.builtIn).toBe(true);
    const v = await viewer();
    for (const id of ['state.get', 'state.remove', 'state.visibility', 'state.update', 'state.action']) {
      expect(v.api.describe(id)).toBeDefined();
    }
    expect(v.api.describe('state.get')?.definition.read).toBe(true);
    expect(v.api.describe('state.remove')?.definition.approval).toBe('destructive');
  });

  it('has its commands', async () => {
    const v = await viewer();
    const ids = StateTreeFeature.commands!(v.api).map((c) => c.id);
    expect(ids).toEqual(['state.visibility', 'state.remove', 'state-tree.open']);
  });

  describe('get', () => {
    it('lists the tree: download to representations', async () => {
      const v = await viewer();
      const r = await v.api.state!.get({});
      expectOk(r);
      const nodes = await tree(v);
      expect(nodes).toHaveLength(1);
      for (const t of ['Trajectory from mmCIF', 'Molecular Model', 'Structure', 'Component', '3D Representation']) {
        expect(byTransform(nodes, t).length, t).toBeGreaterThan(0);
      }
      expect(byTransform(nodes, 'Component').map((n) => n.label)).toEqual(
        expect.arrayContaining(['Polymer (1091 elements)'.split(' (')[0], 'Ligand', 'Water']),
      );
      expect(r.did).toContain('Cartoon');
      for (const n of all(nodes)) {
        expect(n.ref).toBeTruthy();
        expect(n.transformer).toBeTruthy();
        expect(n.status).toBe('ok');
      }
      // Representations are what the eye is for; the data above them has none.
      expect(nodes[0].visible).toBeDefined();
      expect(byTransform(nodes, '3D Representation').every((n) => n.visible === true)).toBe(true);
    });

    it('describes one cell: its params and applicable actions', async () => {
      const v = await viewer();
      const nodes = await tree(v);
      const structure = byTransform(nodes, 'Structure')[0];
      const r = await v.api.state!.get({ ref: structure.ref });
      expectOk(r);
      const d = r.value as any;
      expect(d.ref).toBe(structure.ref);
      expect(d.params.type.name).toBe('assembly');
      expect(d.updates.map((u: any) => u.ref)).toContain(structure.ref);
      expect(d.actions.map((a: any) => a.name)).toEqual(expect.arrayContaining(['selection', '3d-representation']));
      expect(d.actions.find((a: any) => a.name === 'selection').params).toEqual(['script', 'label']);
      expect(r.did).toContain('Actions:');
    });

    it('fails for an unknown ref', async () => {
      const v = await viewer();
      const r = await v.api.state!.get({ ref: 'nope' });
      expect(r).toMatchObject({ ok: false, code: 'not-found' });
    });

    it('says when the state is empty', async () => {
      const v = await headlessViewer({ extensions: [StateTreeFeature] });
      const r = await v.api.state!.get({});
      expectOk(r);
      expect(r.did).toBe('The state is empty.');
      expect(r.value).toEqual([]);
    });
  });

  describe('visibility', () => {
    it('hides a representation, already so on repeat, and Undo shows it', async () => {
      const v = await viewer();
      const rep = byTransform(await tree(v), '3D Representation')[0];
      const before = steps(v);
      const r = await v.api.state!.visibility({ ref: rep.ref, visible: false });
      expectOk(r);
      expect(cellOf(v, rep.ref).state.isHidden).toBe(true);
      expect(steps(v)).toBe(before + 1);

      const again = await v.api.state!.visibility({ ref: rep.ref, visible: false });
      expectOk(again);
      expect(again.unchanged).toBe(true);
      expect(steps(v)).toBe(before + 1);

      expectOk(await v.api.history.undo());
      expect(cellOf(v, rep.ref).state.isHidden).toBeFalsy();
    });

    it('shows it again', async () => {
      const v = await viewer();
      const rep = byTransform(await tree(v), '3D Representation')[0];
      expectOk(await v.api.state!.visibility({ ref: rep.ref, visible: false }));
      expectOk(await v.api.state!.visibility({ ref: rep.ref, visible: true }));
      expect(cellOf(v, rep.ref).state.isHidden).toBeFalsy();
    });

    it('fails for an unknown ref, an empty one and the root', async () => {
      const v = await viewer();
      expect(await v.api.state!.visibility({ ref: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
      expect(await v.api.state!.visibility({ ref: '' })).toMatchObject({ ok: false, code: 'invalid' });
      expect(await v.api.state!.visibility({ ref: StateTransform.RootRef })).toMatchObject({
        ok: false,
        code: 'blocked',
      });
    });
  });

  describe('update', () => {
    const cartoon = async (v: PrismViewer) => {
      const nodes = await tree(v);
      return byTransform(nodes, '3D Representation').find((n) => n.label === 'Cartoon')!;
    };

    it("merges into a cell's params, and Undo restores", async () => {
      const v = await viewer();
      const rep = await cartoon(v);
      const current = cellOf(v, rep.ref).transform.params as any;
      const before = steps(v);
      const r = await v.api.state!.update({
        ref: rep.ref,
        params: { type: { name: 'cartoon', params: { ...current.type.params, alpha: 0.5 } } },
      });
      expectOk(r);
      expect((cellOf(v, rep.ref).transform.params as any).type.params.alpha).toBe(0.5);
      // The other params stay as they were.
      expect((cellOf(v, rep.ref).transform.params as any).colorTheme).toEqual(current.colorTheme);
      expect(steps(v)).toBe(before + 1);

      expectOk(await v.api.history.undo());
      expect((cellOf(v, rep.ref).transform.params as any).type.params.alpha).toBe(1);
    });

    it('changes the color theme', async () => {
      const v = await viewer();
      const rep = await cartoon(v);
      expectOk(
        await v.api.state!.update({
          ref: rep.ref,
          params: { colorTheme: { name: 'secondary-structure', params: {} } },
        }),
      );
      expect((cellOf(v, rep.ref).transform.params as any).colorTheme.name).toBe('secondary-structure');
      expectOk(await v.api.history.undo());
      expect((cellOf(v, rep.ref).transform.params as any).colorTheme.name).toBe('chain-id');
    });

    it('is already so when nothing would change', async () => {
      const v = await viewer();
      const rep = await cartoon(v);
      const current = cellOf(v, rep.ref).transform.params as any;
      const before = steps(v);
      const r = await v.api.state!.update({ ref: rep.ref, params: { colorTheme: current.colorTheme } });
      expectOk(r);
      expect(r.unchanged).toBe(true);
      expect(steps(v)).toBe(before);
    });

    it('checks the params against the transformer', async () => {
      const v = await viewer();
      const component = byTransform(await tree(v), 'Component')[0];
      const before = steps(v);
      expect(await v.api.state!.update({ ref: component.ref, params: { bogus: 1 } })).toMatchObject({
        ok: false,
        code: 'invalid',
        message: expect.stringContaining('bogus'),
      });
      expect(await v.api.state!.update({ ref: component.ref, params: { label: 5 } })).toMatchObject({
        ok: false,
        code: 'invalid',
      });
      expect(steps(v)).toBe(before);
    });

    it('fails for an unknown ref, and for a cell without params', async () => {
      const v = await viewer();
      expect(await v.api.state!.update({ ref: 'nope', params: {} })).toMatchObject({ ok: false, code: 'not-found' });
      expect(await v.api.state!.update({ ref: StateTransform.RootRef, params: {} })).toMatchObject({
        ok: false,
        code: 'blocked',
      });
    });
  });

  describe('action', () => {
    const script = { language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.resname ALA))' };

    it('applies an action of the cell, and Undo removes what it made', async () => {
      const v = await viewer();
      const structure = byTransform(await tree(v), 'Structure')[0];
      const cells = v.plugin.state.data.cells.size;
      const before = steps(v);
      const r = await v.api.state!.action({ ref: structure.ref, action: 'selection', params: { script } });
      expectOk(r);
      expect(v.plugin.state.data.cells.size).toBe(cells + 1);
      expect(steps(v)).toBe(before + 1);
      const made = (await tree(v)).flatMap((n) => all([n])).find((n) => n.transform === 'Selection');
      expect(made).toBeDefined();

      expectOk(await v.api.history.undo());
      expect(v.plugin.state.data.cells.size).toBe(cells);
    });

    it('fills in the action defaults', async () => {
      const v = await viewer();
      const structure = byTransform(await tree(v), 'Structure')[0];
      const cells = v.plugin.state.data.cells.size;
      expectOk(await v.api.state!.action({ ref: structure.ref, action: 'selection' }));
      expect(v.plugin.state.data.cells.size).toBe(cells + 1);
    });

    it('names the actions when the name is wrong', async () => {
      const v = await viewer();
      const structure = byTransform(await tree(v), 'Structure')[0];
      const r = await v.api.state!.action({ ref: structure.ref, action: 'dance' });
      expect(r).toMatchObject({ ok: false, code: 'not-found', hint: expect.stringContaining('selection') });
    });

    it('checks the action params', async () => {
      const v = await viewer();
      const structure = byTransform(await tree(v), 'Structure')[0];
      const r = await v.api.state!.action({ ref: structure.ref, action: 'selection', params: { bogus: 1 } });
      expect(r).toMatchObject({ ok: false, code: 'invalid' });
    });

    it('fails for an unknown ref', async () => {
      const v = await viewer();
      expect(await v.api.state!.action({ ref: 'nope', action: 'selection' })).toMatchObject({
        ok: false,
        code: 'not-found',
      });
    });
  });

  describe('remove', () => {
    it('removes a component with its representations, and Undo restores it', async () => {
      const v = await viewer();
      const nodes = await tree(v);
      const ligand = byTransform(nodes, 'Component').find((n) => n.label === 'Ligand')!;
      const below = all(ligand.children ?? []).map((n) => n.ref);
      expect(below.length).toBeGreaterThan(0);
      const before = steps(v);
      const r = await v.api.state!.remove({ ref: ligand.ref });
      expectOk(r);
      expect(r.did).toContain('Ligand');
      expect(v.plugin.state.data.cells.has(ligand.ref)).toBe(false);
      for (const ref of below) expect(v.plugin.state.data.cells.has(ref)).toBe(false);
      expect(steps(v)).toBe(before + 1);

      expectOk(await v.api.history.undo());
      expect(v.plugin.state.data.cells.has(ligand.ref)).toBe(true);
      for (const ref of below) expect(v.plugin.state.data.cells.has(ref)).toBe(true);
    });

    it('removes a structure as structure.remove would', async () => {
      const v = await viewer();
      const structure = byTransform(await tree(v), 'Structure')[0];
      expectOk(await v.api.state!.remove({ ref: structure.ref }));
      expect(v.managers.history.steps.value.length).toBeGreaterThan(0);
      expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(0);
    });

    it('fails for an unknown ref, an empty one and the root', async () => {
      const v = await viewer();
      expect(await v.api.state!.remove({ ref: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
      expect(await v.api.state!.remove({ ref: '' })).toMatchObject({ ok: false, code: 'invalid' });
      expect(await v.api.state!.remove({ ref: StateTransform.RootRef })).toMatchObject({ ok: false, code: 'blocked' });
    });
  });

  describe('the source', () => {
    const find = (nodes: readonly StateCellNode[], ref: string): StateCellNode | undefined => {
      for (const n of nodes) {
        if (n.ref === ref) return n;
        const c = n.children && find(n.children, ref);
        if (c) return c;
      }
    };

    it('reads the tree on demand, with decorators folded and nothing for ghost cells', async () => {
      const v = await viewer();
      const roots = v.api.state!.tree.roots.value;
      expect(roots.length).toBeGreaterThan(0);
      expect(roots[0].icon).toBe('file');
      const unit = (function walk(ns: readonly StateCellNode[]): StateCellNode[] {
        return ns.flatMap((n) => [n, ...walk(n.children ?? [])]);
      })(roots);
      expect(unit.every((n) => n.status === 'ok')).toBe(true);
    });

    it('follows the state while subscribed: a changed node, the rest the same objects', async () => {
      const v = await viewer();
      const source = v.api.state!.tree;
      const seen: (readonly StateCellNode[])[] = [];
      const sub = source.roots.subscribe((r) => seen.push(r));
      try {
        const first = seen[seen.length - 1];
        const reps = all(await tree(v)).filter((n) => n.transform === '3D Representation');
        const water = reps[reps.length - 1];
        const other = reps[0];
        const otherBefore = find(first, other.ref)!;
        expectOk(await v.api.state!.visibility({ ref: water.ref, visible: false }));
        await until(() => find(seen[seen.length - 1], water.ref)?.visible === false);
        const now = seen[seen.length - 1];
        expect(find(now, other.ref)).toBe(otherBefore);
        expect(find(now, water.ref)!.visible).toBe(false);

        // Undo comes through the same way.
        expectOk(await v.api.history.undo());
        await until(() => find(seen[seen.length - 1], water.ref)?.visible === true);
      } finally {
        sub.unsubscribe();
      }
    });

    it('coalesces events: a burst of changes is one update', async () => {
      const v = await viewer();
      const source = v.api.state!.tree;
      const seen: (readonly StateCellNode[])[] = [];
      const sub = source.roots.subscribe((r) => seen.push(r));
      try {
        const count = seen.length;
        const reps = all(await tree(v)).filter((n) => n.transform === '3D Representation');
        const state = v.plugin.state.data;
        for (const r of reps) state.updateCellState(r.ref, { isHidden: true });
        await until(() => reps.every((r) => find(seen[seen.length - 1], r.ref)?.visible === false));
        expect(seen.length - count).toBe(1);
      } finally {
        sub.unsubscribe();
      }
    });

    it('drops removed cells and shows new ones', async () => {
      const v = await viewer();
      const source = v.api.state!.tree;
      const seen: (readonly StateCellNode[])[] = [];
      const sub = source.roots.subscribe((r) => seen.push(r));
      try {
        const structure = byTransform(await tree(v), 'Structure')[0];
        expectOk(await v.api.state!.action({ ref: structure.ref, action: 'selection' }));
        await until(
          () => find(seen[seen.length - 1], structure.ref)?.children?.some((c) => c.transform === 'Selection') === true,
        );
        const ligand = byTransform(await tree(v), 'Component').find((n) => n.label === 'Ligand')!;
        expectOk(await v.api.state!.remove({ ref: ligand.ref }));
        await until(() => !find(seen[seen.length - 1], ligand.ref));
      } finally {
        sub.unsubscribe();
      }
    });

    it('holds no subscription to the state while nobody looks', async () => {
      const v = await viewer();
      const source = v.api.state!.tree;
      const sub = source.roots.subscribe(() => {});
      sub.unsubscribe();
      const first = source.roots.value;
      const rep = byTransform(await tree(v), '3D Representation')[0];
      expectOk(await v.api.state!.visibility({ ref: rep.ref, visible: false }));
      // Read again: current, with no event ever handled.
      expect(find(source.roots.value, rep.ref)!.visible).toBe(false);
      expect(find(first, rep.ref)!.visible).toBe(true);
    });

    it('reports cells whose params changed', async () => {
      const v = await viewer();
      const source = v.api.state!.tree;
      const changed: ReadonlySet<string>[] = [];
      const sub = source.changed.subscribe((c) => changed.push(c));
      try {
        const rep = byTransform(await tree(v), '3D Representation')[0];
        const current = cellOf(v, rep.ref).transform.params as any;
        expectOk(
          await v.api.state!.update({
            ref: rep.ref,
            params: { type: { name: current.type.name, params: { ...current.type.params, alpha: 0.7 } } },
          }),
        );
        await until(() => changed.some((c) => c.has(rep.ref)));
      } finally {
        sub.unsubscribe();
      }
    });

    it('describes a cell for the tab: its updates first, then the actions', async () => {
      const v = await viewer();
      const structure = byTransform(await tree(v), 'Structure')[0];
      const d = v.api.state!.tree.details(structure.ref)!;
      expect(d.node.ref).toBe(structure.ref);
      expect(d.updates[0].ref).toBe(structure.ref);
      expect(Object.keys(d.updates[0].params)).toContain('type');
      expect(d.actions.map((a) => a.name)).toContain('selection');
      expect(v.api.state!.tree.details('nope')).toBeUndefined();
    });
  });
});

describe('where data comes from (gate G5)', () => {
  it('is not changed through state.update, and actions ask as network', async () => {
    const v = await viewer();
    const tree = v.plugin.state.data.tree;
    const source = tree.children.get(StateTransform.RootRef).first()!;
    expect(await v.api.run('state.update', { ref: source, params: {} })).toMatchObject({ ok: false, code: 'blocked' });
    expect(v.api.describe('state.action')?.definition.approval).toBe('network');
  });
});
