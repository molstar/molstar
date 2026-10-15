/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, crambin, headlessViewer } from '../../_test/harness.js';
import { mapFile, particlesFile } from '../../api/_test/volume-fixtures.js';
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

describe('the scene: volumes and particles', () => {
  it('shows a volume as a top-level node with its representations, next to structures', async () => {
    const v = await headlessViewer({ files: [cbs(), mapFile()] });
    const roots = v.views.scene.roots.value;
    expect(kinds(roots)).toEqual(['structure:1CBS', 'volume:test.ccp4']);
    const volume = roots[1];
    expect(volume.icon).toBe('volume');
    expect(volume.description.value).toMatch(/12 × 12 × 12 · 1\.5 Å/);
    const reps = volume.children.value;
    expect(reps.map((r) => r.kind)).toEqual(['volume-representation']);
    expect(reps[0].label.value).toBe('Isosurface');
    expect(reps[0].description.value).toMatch(/2\.00 σ/);
  });

  it('keeps volume identity across rebuilds and follows an iso value change', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const volume = v.views.scene.roots.value[0];
    const rep = volume.children.value[0];
    await v.api.volume.isovalue({ value: 1, kind: 'relative' });
    expect(v.views.scene.roots.value[0]).toBe(volume);
    expect(volume.children.value[0]).toBe(rep);
    expect(rep.description.value).toMatch(/1\.00 σ/);
    await v.api.load.files({ files: [mapFile('other.ccp4')] });
    expect(v.views.scene.roots.value).toHaveLength(2);
    expect(v.views.scene.roots.value[0]).toBe(volume);
  });

  it('inspects a volume with its grid facts and a volume representation with its form', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const volume = v.views.scene.roots.value[0];
    const facts = new Map(v.views.scene.inspect(volume).facts ?? []);
    expect(facts.get('Grid')).toBe('12 × 12 × 12');
    expect(facts.get('Voxel')).toBe('1.5 Å');
    expect(facts.get('Mean · σ')).toMatch(/ · /);
    expect(facts.get('Range')).toMatch(/ … /);
    expect(facts.get('Source')).toBeDefined();

    const rep = volume.children.value[0];
    const inspection = v.views.scene.inspect(rep);
    expect(inspection.form).toBeDefined();
    expect(inspection.slot).toBe('inspector:volume-representation');
  });

  it('applies a volume representation form through volume.update, as the person', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const rep = v.views.scene.roots.value[0].children.value[0];
    const form = v.views.scene.inspect(rep).form!;
    const current = form.values.value as any;
    form.set('type', { ...current.type, params: { ...current.type.params, alpha: 0.25 } } as never);
    await new Promise((r) => setTimeout(r, 400));
    const cell = v.plugin.managers.volume.hierarchy.current.volumes[0].representations[0].cell;
    expect((cell.transform.params as any).type.params.alpha).toBe(0.25);
    const step = v.managers.history.steps.value[0]!;
    expect(step.origin.kind).toBe('person');
  });

  it('acts on a volume through the API: the eye and the remove are history steps', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const volume = v.views.scene.roots.value[0];
    const rep = volume.children.value[0];

    rep.toggleVisibility();
    await new Promise((r) => setTimeout(r, 20));
    expect(rep.isHidden.value).toBe(true);
    expect(volume.isHidden.value).toBe(false);
    expect(v.managers.history.steps.value[0]?.origin.kind).toBe('person');

    volume.toggleVisibility();
    await new Promise((r) => setTimeout(r, 20));
    expect(volume.isHidden.value).toBe(true);

    v.views.scene.pick(rep);
    volume.remove();
    await new Promise((r) => setTimeout(r, 50));
    expect(v.views.scene.roots.value).toHaveLength(0);
    expect(v.views.scene.picked.value).toBeUndefined();

    await v.api.history.undo();
    expect(v.views.scene.roots.value.map((n) => n.label.value)).toEqual(['test.ccp4']);
  });

  it('removes one volume representation from its row', async () => {
    const v = await headlessViewer({ files: [mapFile()] });
    const volume = v.views.scene.roots.value[0];
    volume.children.value[0].remove();
    await new Promise((r) => setTimeout(r, 50));
    expect(volume.children.value).toHaveLength(0);
    expect(v.views.scene.roots.value).toHaveLength(1);
  });

  it('shows particle lists with their representations, facts and form', async () => {
    const v = await headlessViewer({ files: [particlesFile()] });
    const list = v.views.scene.roots.value[0];
    expect(list.kind).toBe('particles');
    expect(list.description.value).toBe('5 particles');
    const rep = list.children.value[0];
    expect(rep.kind).toBe('particles-representation');
    expect(new Map(v.views.scene.inspect(list).facts ?? []).get('Particles')).toBe('5');
    expect(v.views.scene.inspect(rep).form).toBeDefined();

    rep.toggleVisibility();
    await new Promise((r) => setTimeout(r, 20));
    expect(rep.isHidden.value).toBe(true);

    list.remove();
    await new Promise((r) => setTimeout(r, 50));
    expect(v.views.scene.roots.value).toHaveLength(0);
  });
});
