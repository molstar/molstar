/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree } from '@molstar/core/state';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer } from '../../_test/harness.js';
import { checkpointSize, createCheckpointPort } from '../checkpoint.js';

// The history's checkpoints (sequence.md §7.2, §7.4): what they hold, what they measure, and what Undo and Redo give
// back from them.

/** The full size the budget estimates: the checkpoint's JSON without data payloads. */
const fullSize = (value: unknown) =>
  JSON.stringify(value, (_k, v) => (typeof v === 'string' && v.length > 1024 ? '' : ArrayBuffer.isView(v) ? '' : v))
    .length;

/** The data tree as JSON: what a restore must give back. */
const treeOf = (v: Awaited<ReturnType<typeof headlessViewer>>) => JSON.stringify(StateTree.toJSON(v.plugin.state.data.tree));

describe('checkpoints', () => {
  it('hold the data tree as the plugin serializes it, and follow each change', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const port = createCheckpointPort(v.plugin);
    const a = port.checkpoint();
    expect(a.plugin.data?.tree).toEqual(StateTree.toJSON(v.plugin.state.data.tree));
    expect(a.plugin.camera).toBeUndefined();
    expect(a.plugin.structureSelection).toBeUndefined();

    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'hydrophobicity' }));
    const b = port.checkpoint();
    expect(b.plugin.data?.tree).toEqual(StateTree.toJSON(v.plugin.state.data.tree));
    expect(b.plugin.data?.tree).not.toEqual(a.plugin.data?.tree);
  });

  it('share the transforms a step did not change, and see state flags changed in place', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const port = createCheckpointPort(v.plugin);
    const a = port.checkpoint();
    expectOk(await v.api.representation.color({ target: 'polymer', theme: 'hydrophobicity' }));
    const b = port.checkpoint();
    const byRef = new Map(a.plugin.data!.tree.transforms.map((t) => [t.ref, t]));
    const shared = b.plugin.data!.tree.transforms.filter((t) => byRef.get(t.ref) === t);
    const changed = b.plugin.data!.tree.transforms.filter((t) => byRef.get(t.ref) !== t);
    expect(changed.length).toBeGreaterThan(0);
    expect(shared.length).toBeGreaterThan(changed.length);
    for (const t of changed) expect(t.version).not.toBe(byRef.get(t.ref)?.version);

    // The state's transient tree merges a second flag update into the same transform object.
    const ref = [...v.plugin.state.data.cells.values()].find((c) =>
      c.transform.transformer.id.endsWith('structure-representation-3d'),
    )!.transform.ref;
    const hidden = (c: ReturnType<typeof port.checkpoint>) =>
      !!c.plugin.data!.tree.transforms.find((t) => t.ref === ref)?.state?.isHidden;
    v.plugin.state.data.updateCellState(ref, { isHidden: true });
    expect(hidden(port.checkpoint())).toBe(true);
    v.plugin.state.data.updateCellState(ref, { isHidden: false });
    expect(hidden(port.checkpoint())).toBe(false);
    v.plugin.state.data.updateCellState(ref, { isHidden: true });
    const c = port.checkpoint();
    expect(hidden(c)).toBe(true);
    expect(c.plugin.data?.tree).toEqual(StateTree.toJSON(v.plugin.state.data.tree));
  });

  it('measure the JSON size without payloads', async () => {
    const v = await headlessViewer({ files: [cbs(), crambin()] });
    const port = createCheckpointPort(v.plugin);
    const c = port.checkpoint();
    const size = port.size!(c);
    expect(size).toBeGreaterThan(1000);
    expect(Math.abs(size - fullSize(c.plugin)) / fullSize(c.plugin)).toBeLessThan(0.02);
    expect(Math.abs(checkpointSize(c) - size)).toBeLessThan(0.02 * size);
    // A change is measured: another style changes the size.
    expectOk(await v.api.styles.apply({ style: 'spacefill' }));
    const d = port.checkpoint();
    expect(port.size!(d)).not.toBe(size);
    expect(Math.abs(port.size!(d) - fullSize(d.plugin)) / fullSize(d.plugin)).toBeLessThan(0.02);
  });

  it('undo and redo walk back and forth through several steps', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const states = [treeOf(v)];
    for (const theme of ['hydrophobicity', 'secondary-structure', 'element-symbol'] as const) {
      expectOk(await v.api.representation.color({ target: 'polymer', theme }));
      states.push(treeOf(v));
    }
    expect(new Set(states).size).toBe(4);
    for (let i = 3; i > 0; i--) {
      expectOk(await v.api.history.undo({}));
      expect(treeOf(v)).toEqual(states[i - 1]);
    }
    for (let i = 1; i <= 3; i++) {
      expectOk(await v.api.history.redo({}));
      expect(treeOf(v)).toEqual(states[i]);
    }
    // Undo after redo still restores the step before.
    expectOk(await v.api.history.undo({}));
    expect(treeOf(v)).toEqual(states[2]);
  });
});
