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
import { shownInOrder } from '@molstar/prism-ui/testing/virtual-list';
import { describe, expect, it } from 'vitest';
import { cbs, crambin, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import type { SceneModel } from '../../../model/scene.js';
import { SceneTree } from '../scene-tree.js';
import { componentsOf } from '../../../testing/scene-groups.js';

// The scene tree over the real scene source (its change stream, sources/scene.ts): after each change a person or the
// API makes, its rows are those a tree mounted now makes from scratch.

(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= class {
  observe() {}
  unobserve() {}
  disconnect() {}
};

function mount(scene: SceneModel) {
  const box = document.createElement('div');
  document.body.append(box);
  const dispose = render(() => <SceneTree scene={scene} />, box);
  flush();
  const rows = () => {
    flush();
    return shownInOrder(box, '[role=treeitem]').map((e) =>
      [
        e.dataset.treeKey,
        e.getAttribute('aria-level'),
        e.getAttribute('aria-posinset'),
        e.getAttribute('aria-setsize'),
        e.getAttribute('aria-expanded'),
      ].join(' '),
    );
  };
  return {
    rows,
    dispose: () => {
      dispose();
      box.remove();
    },
  };
}

describe('the scene tree over the real scene', () => {
  it('matches a tree made from scratch after each change', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const scene = v.views.scene;
    const t = mount(scene);
    const same = (what: string) => {
      const fresh = mount(scene);
      expect(t.rows(), what).toEqual(fresh.rows());
      fresh.dispose();
    };
    same('loaded');
    const s = scene.roots.value[0];
    const contents = s.children.value[0];
    contents.toggleExpanded();
    same('contents opened');
    const entity = contents.children.value[0];
    entity.toggleExpanded();
    same('an entity opened');
    const component = componentsOf(s).find((c) => c.kind === 'component')!;
    component.toggleExpanded();
    same('a component opened');
    const count = componentsOf(s).length;
    component.remove();
    await until(() => componentsOf(s).length === count - 1);
    same('a component removed');
    expectOk(await v.api.history.undo());
    await until(() => componentsOf(s).length === count);
    same('undone');
    expectOk(await v.api.load.files({ files: [crambin()] }));
    same('another structure');
    s.toggleExpanded();
    same('the structure closed');
    t.dispose();
  });
});
