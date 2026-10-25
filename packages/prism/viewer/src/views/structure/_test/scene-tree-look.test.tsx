// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { flush } from '@molstar/prism-ui/dom/scope';
import { describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import { classes } from '../_generated/scene-tree.module.css.js';
import { SceneTree } from '../scene-tree.js';

// The scene tree's look marks over the real scene (api.md §2, "Look"): the sliders glyph on a row with its own look,
// faint on the rows above, its tooltip saying what; gone when the default comes to equal it.

(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= class {
  observe() {}
  unobserve() {}
  disconnect() {}
};

describe('the look marks in the scene tree', () => {
  it('marks the component with its own material, and the structure above it; follows the defaults', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const box = document.createElement('div');
    document.body.append(box);
    const dispose = render(() => <SceneTree scene={v.views.scene} />, box);
    const marks = () => {
      flush();
      return [...box.querySelectorAll<HTMLElement>(`.${classes.look}`)]
        .filter((m) => m.style.display !== 'none')
        .map((m): [string | undefined, string | undefined, boolean] => [
          m.closest('[data-tree-row]')?.textContent?.split(' ')[0],
          m.dataset.tip,
          m.classList.contains(classes.lookInside),
        ]);
    };
    try {
      await until(() => v.views.scene.roots.value.length > 0);
      expect(marks()).toEqual([]);
      const ligand = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
        (c) => c.cell.obj?.label === 'Ligand',
      )!;
      const Metallic = { metalness: 1, roughness: 0.6, bumpiness: 0 };
      expectOk(await v.api.look.set({ node: ligand.cell.transform.ref, material: Metallic }));
      await until(() => marks().length > 0);
      const shown = marks();
      expect(shown.find(([row]) => row?.startsWith('Ligand'))).toEqual([
        expect.stringMatching(/^Ligand/),
        'A representation here has its own look',
        true,
      ]);
      expect(shown.some(([row]) => row?.startsWith('1CBS'))).toBe(true);

      // the default comes to equal it: no longer its own (an options change that rewrites nothing)
      expectOk(await v.api.look.set({ material: Metallic }));
      await until(() => marks().length === 0);
    } finally {
      dispose();
      box.remove();
    }
  });
});
