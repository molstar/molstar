/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer, until } from '../../_test/harness.js';

// A component's actions (its row menu, its inspector) open their commands in the command menu for it: the component
// and the structure filled in, a preset component's keyword as the target, the rest to choose there.

describe("a component's actions", () => {
  it('open their commands with the component filled in', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    await until(() => v.views.scene.roots.value.length > 0);
    const polymer = v.views.scene.roots.value[0].children.value.find((n) => n.label.value === 'Polymer')!;
    const run = (id: string) => polymer.menu!.find((m) => m.id === id)!.run!();
    const active = () => v.managers.commands.active.value;

    run('add-repr');
    await until(() => !!active());
    expect(active()!.command.id).toBe('representation.set');
    expect(active()!.state.value.values).toMatchObject({
      component: polymer.id,
      structure: v.views.scene.roots.value[0].id,
    });

    run('coloring');
    await until(() => active()?.command.id === 'representation.color-theme');
    expect(active()!.state.value.values).toMatchObject({ target: 'polymer' });

    run('modify');
    await until(() => active()?.command.id === 'component.modify');
    expect(active()!.state.value.values).toMatchObject({ component: polymer.id, target: 'selection' });
  });
});
