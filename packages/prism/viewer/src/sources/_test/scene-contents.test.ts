/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer, until } from '../../_test/harness.js';
import type { SceneNode } from '../../model/scene.js';

// The contents rows (README decision 18): an entity's or a chain's eye hides it and shows it again, the row says which,
// and hovering it highlights its atoms.

describe('the contents rows', () => {
  it('show what is hidden, toggle back, and highlight on hover', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const contents = v.views.scene.roots.value[0].children.value.find((c) => c.kind === 'contents')!;
    contents.load?.();
    const entity = () =>
      contents.children.value.find((c) => c.kind === 'entity' && /RETINOIC/i.test(c.label.value)) as SceneNode;
    expect(entity().isHidden.value).toBe(false);

    entity().toggleVisibility();
    await until(() => entity().isHidden.value);
    expect(entity().children.value.every((c) => c.isHidden.value)).toBe(true);
    const components = () => v.plugin.managers.structure.hierarchy.current.structures[0].components.length;
    expect(components()).toBe(2);

    entity().toggleVisibility();
    await until(() => !entity().isHidden.value);
    await until(() => components() === 3);

    entity().highlight(true);
    await until(() => v.sources.pointer.hoverLines.value.length > 0);
    entity().highlight(false);
    await until(() => v.sources.pointer.hoverLines.value.length === 0);
  });
});
