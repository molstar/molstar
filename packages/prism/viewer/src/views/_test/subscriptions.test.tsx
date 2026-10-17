// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { render } from '@molstar/prism-ui/dom/render';
import { bind, flush, flushFrame } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ReactiveProbe } from '@molstar/prism-ui/reactive/state';
import { makeScene } from '../../_stories/performance-fixtures.js';
import type { StateTreeModel, StateTreeNode } from '../../model/state-tree.js';
import { StateTreeTab } from '../state-tree/state-tree-tab.js';
import { SceneTree } from '../structure/scene-tree.js';

// Subscriptions scale with what is shown, not with the model (view-models.md §3): a tree view's live listeners are a
// constant and some per row on screen, the same for a model four times larger, and within its budget
// (subscription-budget.json, shrink-only: `SUBS_UPDATE=1` writes today's after an intended change). A view that listens
// to every item (the positive control) is caught. The listeners are counted at `listen` (`ReactiveProbe.listened`).

(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= class {
  observe() {}
  unobserve() {}
  disconnect() {}
};

let live = 0;
ReactiveProbe.on = true;
ReactiveProbe.listened = () => {
  live++;
  return 1;
};
ReactiveProbe.unlistened = () => {
  live--;
};

afterEach(() => {
  document.body.innerHTML = '';
});

/** The listeners a view adds when mounted, and the rows it shows. */
function mounted(view: () => Node) {
  const box = document.createElement('div');
  document.body.append(box);
  const before = live;
  const dispose = render(view, box);
  flush();
  return { listeners: live - before, rows: box.querySelectorAll('[role=treeitem]').length, box, dispose };
}

/** 1,000 components open (4,000 rows): `components` components. */
function openScene(components: number) {
  const fx = makeScene(components, 3);
  for (const c of fx.components) c.expanded.set(true);
  return fx;
}

/** `structures` structures of 4 models of 49 cells (201 cells each), all open. */
function stateTree(structures: number): StateTreeModel {
  const roots: StateTreeNode[] = [];
  for (let s = 0; s < structures; s++)
    roots.push({
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
          label: `Component ${c}`,
          icon: 'cartoon',
          transform: 'Component',
          visible: true,
        })),
      })),
    });
  return {
    roots: observableState<readonly StateTreeNode[]>(roots).asReadonly,
    selected: observableState<string | undefined>(undefined).asReadonly,
    select: () => {},
    inspect: () => undefined,
    toggleVisible: () => {},
    remove: () => {},
  };
}

const BudgetFile = path.join(import.meta.dirname, 'subscription-budget.json');
const budget = (): Record<string, number> => JSON.parse(fs.readFileSync(BudgetFile, 'utf8'));
const measured: Record<string, number> = {};

/** The same listeners for a model 4× larger, the rows on screen the same, and within the budget. */
function scales(name: string, small: { listeners: number; rows: number }, large: { listeners: number; rows: number }) {
  expect(large.rows, `${name}: the rows on screen`).toBe(small.rows);
  expect(large.listeners, `${name}: listeners for a model 4× larger`).toBe(small.listeners);
  measured[name] = large.listeners;
  if (process.env.SUBS_UPDATE) return;
  expect(large.listeners, `${name}: over its budget (subscription-budget.json)`).toBeLessThanOrEqual(budget()[name]);
}

describe('subscriptions scale with what is shown', () => {
  it('the scene tree: one listener for the structure (its change stream), the rest per row', () => {
    // Made before counting: the fixture's own listeners (its change stream's) are not the view's.
    const small = openScene(250).scene;
    const large = openScene(1000).scene;
    const a = mounted(() => <SceneTree scene={small} />);
    a.dispose();
    const b = mounted(() => <SceneTree scene={large} />);
    b.dispose();
    scales('scene tree', a, b);
  });

  it('the state tree', () => {
    const small = stateTree(25);
    const large = stateTree(100);
    const a = mounted(() => <StateTreeTab model={small} />);
    a.dispose();
    const b = mounted(() => <StateTreeTab model={large} />);
    b.dispose();
    scales('state tree', a, b);
  });

  it("a thumb's drag tells the scene's nodes nothing (no model state changes while the rows move)", () => {
    const fx = openScene(1000);
    const t = mounted(() => <SceneTree scene={fx.scene} />);
    const nodeStates = new Set<object>();
    for (const n of [fx.structure, ...fx.components]) nodeStates.add(n.expanded).add(n.children).add(n.label);
    const viewport = t.box.querySelector<HTMLElement>('[data-axis] > div')!;
    let top = 0;
    Object.defineProperty(viewport, 'scrollTop', { get: () => top, set: (v: number) => (top = v), configurable: true });
    let told = 0;
    let modelTold = 0;
    ReactiveProbe.told = (state) => {
      told++;
      if (nodeStates.has(state)) modelTold++;
    };
    try {
      for (let i = 0; i < 200; i++) {
        top = Math.round((((i * 37) % 97) / 97) * 4000 * 26);
        viewport.dispatchEvent(new Event('scroll'));
        flushFrame();
        flush();
      }
    } finally {
      ReactiveProbe.told = undefined;
    }
    // The list's own states move (its position, its window); no node's.
    expect(told).toBeGreaterThan(0);
    expect(modelTold).toBe(0);
    t.dispose();
  });

  it('a view that listens to every item is caught (positive control)', () => {
    const everyLabel = (components: number) => {
      const fx = makeScene(components, 3);
      return () => {
        const el = document.createElement('div');
        // Not the narrowest source: a listener per component, shown or not.
        for (const c of fx.components) bind(c.node.label, (l) => (el.title = l));
        return el;
      };
    };
    const a = mounted(everyLabel(250));
    a.dispose();
    const b = mounted(everyLabel(1000));
    b.dispose();
    expect(b.listeners).toBeGreaterThan(a.listeners);
  });

  it('writes the budget (SUBS_UPDATE=1)', () => {
    if (!process.env.SUBS_UPDATE) return;
    fs.writeFileSync(BudgetFile, `${JSON.stringify(measured, null, 2)}\n`);
  });
});
