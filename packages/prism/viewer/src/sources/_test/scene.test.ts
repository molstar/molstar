/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, headlessViewer } from '../../_test/harness.js';
import type { SceneNode } from '../../model/scene.js';

// The scene source (plan 2.2): the Structure tab's tree over a real plugin.

const kinds = (nodes: readonly SceneNode[]) => nodes.map((n) => `${n.kind}:${n.label.value}`);

describe('the scene', () => {
  it('shows each structure with its contents and components, by entry name', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const roots = v.views.scene.roots.value;
    expect(kinds(roots)).toEqual(['structure:1CBS']);
    const children = roots[0].children.value;
    expect(children[0].kind).toBe('contents');
    expect(children.filter((c) => c.kind === 'component').length).toBeGreaterThan(0);
    const component = children.find((c) => c.kind === 'component')!;
    expect(component.children.value[0].kind).toBe('representation');
  });

  it('keeps node identity across rebuilds and follows changes', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const before = v.views.scene.roots.value[0];
    await v.api.load.files({ files: [crambin()] });
    const after = v.views.scene.roots.value;
    expect(after).toHaveLength(2);
    expect(after[0]).toBe(before);
    await v.api.structure.visibility({ structure: '1cbs', visible: false });
    expect(before.isHidden.value).toBe(true);
  });

  it('computes contents when first expanded: entities, then chains', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const contents = v.views.scene.roots.value[0].children.value[0];
    expect(contents.children.value).toHaveLength(0);
    contents.toggleExpanded();
    const entities = contents.children.value;
    expect(entities.length).toBeGreaterThan(1);
    expect(entities[0].children.value[0].kind).toBe('chain');
  });

  it('acts through the API: a row eye is a history step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const s = v.views.scene.roots.value[0];
    s.toggleVisibility();
    await new Promise((r) => setTimeout(r, 10));
    expect(v.managers.history.steps.value[0].label).toMatch(/Hid 1CBS/);
    expect(v.managers.history.steps.value[0].origin.kind).toBe('person');
  });

  it('inspects a structure with facts and a representation with its form', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const s = v.views.scene.roots.value[0];
    const facts = v.views.scene.inspect(s).facts ?? [];
    expect(facts.find(([k]) => k === 'Entry')?.[1]).toBe('1CBS');
    const repr = s.children.value.find((c) => c.kind === 'component')!.children.value[0];
    expect(v.views.scene.inspect(repr).form).toBeDefined();
  });

  it('drops nodes of what is gone', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    v.views.scene.pick(v.views.scene.roots.value[0]);
    await v.api.session.clear();
    expect(v.views.scene.roots.value).toHaveLength(0);
    expect(v.views.scene.picked.value).toBeUndefined();
  });
});
