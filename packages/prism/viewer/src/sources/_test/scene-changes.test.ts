/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, until } from '../../_test/harness.js';
import type { SceneChanges } from '../../model/scene-changes.js';
import { componentsOf } from '../../testing/scene-groups.js';
import { componentsNodeId } from '../../model/scene.js';

// The scene's structure as one change stream (`SceneModel.changes`, model/scene-changes.ts) over a real plugin: the nodes
// whose children or expansion changed, once per Mol* transaction or toggle; nothing for what changes no structure.

async function viewer() {
  const v = await headlessViewer({ files: [cbs()] });
  // The Sequence node's tracks come once the sequence's runtime is in (sequence.md §8.1): a change of its own, before
  // the test's.
  const viewers = (await v.views.sequence!.load()).viewers.get();
  await until(() => viewers.tracks.value.length > 0);
  await until(() =>
    v.views.scene.roots.value[0].children.value.some((c) => c.kind === 'sequence' && c.children.value.length > 0),
  );
  const told: SceneChanges[] = [];
  v.views.scene.changes!.listen?.((c) => told.push(c));
  return { v, told };
}

describe('the scene’s change stream', () => {
  it('a node opened: one change, its id and its lazy children (entities, their chains) loaded in it', async () => {
    const { v, told } = await viewer();
    const contents = v.views.scene.roots.value[0].children.value[0];
    expect(contents.kind).toBe('contents');
    const version = v.views.scene.changes!.value.version;
    contents.toggleExpanded();
    expect(told).toHaveLength(1);
    expect(told[0].nodes.has(contents.id)).toBe(true);
    expect(told[0].version).toBe(version + 1);
    expect(contents.children.value.length).toBeGreaterThan(0);
  });

  it('a transaction that changes no structure (a visibility) tells nothing', async () => {
    const { v, told } = await viewer();
    const s = v.views.scene.roots.value[0];
    expectOk(await v.api.structure.visibility({ structure: '1cbs', visible: false }));
    await until(() => s.isHidden.value);
    expect(told).toEqual([]);
  });

  it("a component removed: its structure's Components group, in one change", async () => {
    const { v, told } = await viewer();
    const s = v.views.scene.roots.value[0];
    const component = componentsOf(s).find((c) => c.kind === 'component')!;
    const count = componentsOf(s).length;
    component.remove();
    await until(() => componentsOf(s).length === count - 1);
    expect(told.length).toBeGreaterThan(0);
    expect(told.some((c) => c.nodes.has(componentsNodeId(s.id)))).toBe(true);
    // Versions count up one a change.
    for (let i = 1; i < told.length; i++) expect(told[i].version).toBe(told[i - 1].version + 1);
  });

  it('reveal opens the path in one change; a new structure is the roots’ (not a change of a node)', async () => {
    const { v, told } = await viewer();
    const s = v.views.scene.roots.value[0];
    const contents = s.children.value[0];
    contents.load?.();
    // Loaded: one change.
    expect(told).toHaveLength(1);
    told.length = 0;
    const entity = contents.children.value[0];
    expect(v.views.scene.reveal!(entity.id)).toBe(true);
    expect(told).toHaveLength(1);
    expect(told[0].nodes.has(contents.id)).toBe(true);
    const roots = v.views.scene.roots.value;
    told.length = 0;
    expectOk(await v.api.load.files({ files: [crambin()] }));
    expect(v.views.scene.roots.value).not.toBe(roots);
  });
});
