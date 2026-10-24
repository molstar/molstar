/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer, until } from '../../_test/harness.js';
import { componentsOf } from '../../testing/scene-groups.js';

// The focus: one row at the top level after the structures (the tree does not shift as the focus comes and goes),
// holding the focus behavior's components while something is focused; named by its structure when there are several.

describe('the focus row', () => {
  it('is there before, during and after a focus, with the target and the surroundings in it', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const focus = () => v.views.scene.roots.value.find((c) => c.kind === 'focus')!;
    // At the top level, right after the structure, and not under it.
    expect(v.views.scene.roots.value.map((n) => n.kind).slice(0, 2)).toEqual(['structure', 'focus']);
    expect(v.views.scene.roots.value[0].children.value.some((c) => c.kind === 'focus')).toBe(false);
    expect(focus().description.value).toBe('none');
    expect(focus().children.value).toEqual([]);
    const components = () => componentsOf(v.views.scene.roots.value[0]).filter((c) => c.kind === 'component');
    const before = components().length;
    expect(before).toBeGreaterThan(0);

    const r = await v.api.focus.set({ target: 'ligand' });
    expect(r.ok).toBe(true);
    await until(() => focus().children.value.length === 2);
    expect(focus().children.value.map((c) => c.label.value)).toEqual(['Target', 'Surroundings (5 Å)']);
    expect(focus().description.value).toMatch(/\d+ atoms/);
    // The focus parts are not among the structure's own components.
    expect(components().length).toBe(before);

    await v.api.focus.clear({});
    await until(() => focus().children.value.length === 0);
    expect(focus().description.value).toBe('none');
  });

  it('hides and shows all its parts from one eye, as one history step', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const focus = () => v.views.scene.roots.value.find((c) => c.kind === 'focus')!;
    await v.api.focus.set({ target: 'ligand' });
    await until(() => focus().children.value.length === 2);
    const steps = () => v.managers.history.steps.value.length;
    const before = steps();

    focus().toggleVisibility();
    await until(() => focus().children.value.every((c) => c.isHidden.value));
    await until(() => steps() > before);
    expect(steps()).toBe(before + 1);

    expect((await v.api.history.undo()).ok).toBe(true);
    await until(() => focus().children.value.every((c) => !c.isHidden.value));
  });
});
