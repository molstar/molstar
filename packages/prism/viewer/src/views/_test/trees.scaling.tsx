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
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { it } from 'vitest';
import { checkLinear, scaling } from '../../../../ui/src/_test/scaling.js';
import { makeScene, type SceneFixture } from '../../_stories/performance-fixtures.js';
import type { StateTreeNode } from '../../model/state-tree.js';
import { StateTreeTab } from '../state-tree/state-tree-tab.js';
import { SceneTree } from '../structure/scene-tree.js';

// The trees' bulk paths grow about linearly (`test:scaling`, prism-ui's _test/scaling.ts): every component of a scene
// opened, then closed, in one task; a structure of N components closed and opened (a reflow of N siblings); the state
// tree's cells all changed at once (the source's flush after a preset over everything).

(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= class {
  observe() {}
  unobserve() {}
  disconnect() {}
};

interface Mounted {
  fx: SceneFixture;
  box: HTMLElement;
  dispose: () => void;
}

function sceneTree(components: number, open: boolean): Mounted {
  const fx = makeScene(components, 3);
  if (open) for (const c of fx.components) c.expanded.set(true);
  const box = document.createElement('div');
  document.body.append(box);
  const dispose = render(() => <SceneTree scene={fx.scene} />, box);
  flush();
  return { fx, box, dispose };
}
const unmount = <M extends { box: HTMLElement; dispose: () => void }>(m: M) => {
  m.dispose();
  m.box.remove();
};

it('the scene tree: every component opened in one task', () => {
  const s = scaling(
    {
      setup: (n) => sceneTree(n, false),
      run: ({ fx }) => {
        for (const c of fx.components) c.expanded.set(true);
        flush();
      },
      teardown: unmount,
    },
    1000,
  );
  checkLinear('scene tree: expand all', s);
});

it('the scene tree: every component closed in one task', () => {
  const s = scaling(
    {
      setup: (n) => sceneTree(n, true),
      run: ({ fx }) => {
        for (const c of fx.components) c.expanded.set(false);
        flush();
      },
      teardown: unmount,
    },
    1000,
  );
  checkLinear('scene tree: collapse all', s);
});

it('the scene tree: a structure of N components closed and opened (a reflow of N siblings)', () => {
  const s = scaling(
    {
      setup: (n) => sceneTree(n, false),
      run: ({ fx }) => {
        fx.structure.expanded.set(false);
        flush();
        fx.structure.expanded.set(true);
        flush();
      },
      teardown: unmount,
    },
    1000,
  );
  checkLinear('scene tree: reflow of N siblings', s);
});

/** `structures` structures of 4 models of 49 cells, labelled with `tag` (a new tag: every cell a new record). */
function cells(structures: number, tag: string): StateTreeNode[] {
  return Array.from({ length: structures }, (_, s) => ({
    ref: `s${s}`,
    label: `Structure ${s}`,
    icon: 'structure',
    transform: 'Structure',
    children: Array.from({ length: 4 }, (_, m) => ({
      ref: `s${s}m${m}`,
      label: `Model ${m}`,
      icon: 'layers',
      transform: 'Model',
      children: Array.from({ length: 49 }, (_, c) => ({
        ref: `s${s}m${m}c${c}`,
        label: `Component ${c} ${tag}`,
        icon: 'cartoon',
        transform: 'Component',
        visible: true,
      })),
    })),
  }));
}

it('the state tree: every cell changed at once', () => {
  const s = scaling(
    {
      setup: (n) => {
        const roots = observableState<readonly StateTreeNode[]>(cells(n, 'a'));
        const next = cells(n, 'b');
        const box = document.createElement('div');
        document.body.append(box);
        const model = {
          roots: roots.asReadonly,
          selected: observableState<string | undefined>(undefined).asReadonly,
          select: () => {},
          inspect: () => undefined,
          toggleVisible: () => {},
          remove: () => {},
        };
        const dispose = render(() => <StateTreeTab model={model} />, box);
        flush();
        return { roots, next, box, dispose };
      },
      run: ({ roots, next }) => {
        roots.set(next);
        flush();
      },
      teardown: unmount,
    },
    10,
  );
  checkLinear('state tree: every cell changed', s);
});
