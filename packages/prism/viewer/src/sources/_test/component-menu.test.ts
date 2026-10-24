/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import { componentsOf } from '../../testing/scene-groups.js';

// A component's actions (its row menu, its inspector) open their commands in the command menu for it: the component
// and the structure filled in, a preset component's keyword as the target, the rest to choose there; Modify by
// selection acts at once.

describe("a component's actions", () => {
  it('open their commands with the component filled in', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    await until(() => v.views.scene.roots.value.length > 0);
    const polymer = componentsOf(v.views.scene.roots.value[0]).find((n) => n.label.value === 'Polymer')!;
    const run = (id: string) => polymer.menu!.find((m) => m.id === id)!.run!();
    const active = () => v.managers.commands.active.value;

    run('add-repr');
    await until(() => !!active());
    expect(active()!.command.id).toBe('representation.set');
    // Named as the person reads them.
    expect(active()!.state.value.values).toMatchObject({ component: 'Polymer', structure: '1CBS' });

    run('coloring');
    await until(() => active()?.command.id === 'representation.color-theme');
    expect(active()!.state.value.values).toMatchObject({ target: 'polymer' });

    // Modify by selection acts at once (three icons on its row), with no form; unavailable with nothing selected.
    const modify = polymer.menu!.find((m) => m.id === 'modify')!;
    expect(modify.inline).toBe(true);
    expect(modify.items?.map((m) => m.id)).toEqual(['add', 'remove', 'intersect']);
    const disabled = () => {
      let d: boolean | undefined;
      modify.items![0].disabled!.subscribe((x) => (d = x)).unsubscribe();
      return d;
    };
    expect(disabled()).toBe(true);
    expectOk(await v.api.select.set({ target: 'ligand' }));
    expect(disabled()).toBe(false);
  });
});

describe('a component in the inspector', () => {
  it('is renamed from its title, and revealed in the tree by id once made', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    await until(() => v.views.scene.roots.value.length > 0);
    const scene = v.views.scene;
    const polymer = componentsOf(scene.roots.value[0]).find((n) => n.label.value === 'Polymer')!;
    // A preset component (a static one: Mol*'s update used to keep its old label).
    scene.inspect(polymer).rename!('The protein');
    await until(() => polymer.label.value === 'The protein');

    // Made with no representation (v5's "< Create Later >"), then revealed: its structure expanded, it picked.
    expectOk(await v.api.select.set({ target: 'ligand' }));
    const made = await v.api.component.add({ target: 'selection', label: 'Mine', representation: 'none' });
    expectOk(made);
    const ref = (made.value as { ref: string }[])[0].ref;
    await until(() => !!scene.reveal?.(ref));
    expect(scene.picked.value?.label.value).toBe('Mine');
    expect(scene.picked.value?.children.value).toEqual([]);
    expect(scene.roots.value[0].isExpanded.value).toBe(true);
  });
});
