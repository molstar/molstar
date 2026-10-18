/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Show } from '@molstar/prism-ui/dom/control-flow';
import { render } from '@molstar/prism-ui/dom/render';
import { flush, flushFrame, Scope } from '@molstar/prism-ui/dom/scope';
import { domView } from '@molstar/prism-ui/dom/view';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { shownInOrder } from '@molstar/prism-ui/testing/virtual-list';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { SceneModel, SceneNode, SceneNodeKind } from '../../../model/scene.js';
import type { SceneChanges } from '../../../model/scene-changes.js';
import type { ViewerUIContext } from '../../../spec.js';
import { extensionUIPort } from '../../../ui/extensions.js';
import { mockScene } from '../../../testing/scene.js';
import { SceneRow, SceneRowPickEvent, type SceneRowPick } from '../scene-row.js';
import { SceneTree } from '../scene-tree.js';
import { StructureTab } from '../structure-tab.js';

// Fixes of the review of the trees: an extension row picked stays picked while it scrolls; rows made by a reflow are the
// rows of the next flattening; a search follows the scene; an extension row's ARIA when its view renders late; the
// search's debounce pauses with the tab.

afterEach(() => {
  document.body.innerHTML = '';
  vi.useRealTimers();
});

function node(id: string, kind: SceneNodeKind, label: string, children: SceneNode[] = [], expanded = false) {
  const isExpanded = observableState(expanded);
  const kids = observableState<readonly SceneNode[]>(children);
  return {
    id,
    kind,
    icon: 'structure',
    label: observableState(label).asReadonly,
    description: observableState<string | undefined>(undefined).asReadonly,
    isHidden: observableState(false).asReadonly,
    isExpanded: isExpanded.asReadonly,
    children: kids.asReadonly,
    setChildren: (c: SceneNode[]) => kids.set(c),
    toggleVisibility: () => {},
    toggleExpanded: () => isExpanded.set(!isExpanded.value),
    highlight: () => {},
    focus: () => {},
    remove: () => {},
  } as unknown as SceneNode & { setChildren(c: SceneNode[]): void };
}

function mount(options: { count?: number; filter?: string } = {}) {
  const components = Array.from({ length: options.count ?? 300 }, (_, i) =>
    node(`c${i}`, 'component', `Component ${i}`),
  );
  const structure = node('s', 'structure', '1ABC', components, true);
  const changes = observableState<SceneChanges>({ version: 0, nodes: new Set() });
  const scene = {
    roots: observableState([structure]).asReadonly,
    picked: observableState<SceneNode | undefined>(undefined).asReadonly,
    pick: () => {},
    changes: changes.asReadonly,
  } as unknown as SceneModel;
  const change = (...ids: string[]) => changes.set({ version: changes.value.version + 1, nodes: new Set(ids) });
  // A change the tree missed one of (resumed): all of it flattened again.
  const skip = () => changes.set({ version: changes.value.version + 2, nodes: new Set() });
  const shell = new ShellModel();
  const port = extensionUIPort(shell, {} as never, {} as ViewerUIContext);
  const filter = observableState(options.filter ?? '');
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <SceneTree scene={scene} shell={shell} filter={filter.asReadonly} />, container);
  const tree = container.querySelector<HTMLElement>('[role=tree]')!;
  let viewport: HTMLElement | undefined;
  for (let e: HTMLElement | null = tree; e; e = e.parentElement)
    if (e.getAttribute('class')?.includes('viewport')) viewport = e;
  const scrollTo = (index: number) => {
    viewport!.scrollTop = index * 26;
    viewport!.dispatchEvent(new Event('scroll'));
    flushFrame();
  };
  const labels = () => {
    flush();
    return shownInOrder(container, '[data-tree-row]').map((r) => r.querySelector('[id]')?.textContent);
  };
  return { structure, components, filter, change, skip, port, container, tree, scrollTo, labels, dispose };
}

describe('an extension row picked in the tree', () => {
  it('stays picked while it scrolls out of the window, and is gone once its row leaves the tree', () => {
    const t = mount();
    const unregister = t.port.register('x', {
      rows: [
        {
          id: 'x.row.membrane',
          when: (n) => n.kind === 'structure',
          render: domView(() => <SceneRow label="Membrane" body={() => <p>Planes</p>} />),
        },
      ],
    });
    flush();
    let pick: SceneRowPick | undefined;
    t.container.addEventListener(SceneRowPickEvent, (e) => {
      pick = (e as CustomEvent<SceneRowPick>).detail;
      e.preventDefault();
    });
    // The extension's row follows the 300 components: scroll to the end to reach it.
    t.scrollTo(300);
    const host = () => t.container.querySelector<HTMLElement>('[data-scene-row="x.row.membrane"]');
    expect(host()).not.toBeNull();
    host()!.querySelector<HTMLElement>('[data-tree-row]')!.click();
    expect(pick).toBeDefined();
    const gone = vi.fn();
    pick!.gone.listen(gone);

    t.scrollTo(0);
    expect(host()).toBeNull();
    expect(gone).not.toHaveBeenCalled();
    // Back in view: made anew, still the picked one.
    t.scrollTo(300);
    expect(host()).not.toBeNull();
    expect(gone).not.toHaveBeenCalled();

    unregister();
    flush();
    expect(gone).toHaveBeenCalledOnce();
    t.dispose();
  });

  it('is the same row after a reflow by changes and a full flattening', () => {
    const t = mount({ count: 5 });
    t.port.register('x', {
      rows: [
        {
          id: 'x.row.membrane',
          when: (n) => n.kind === 'structure',
          render: domView(() => <SceneRow label="Membrane" />),
        },
      ],
    });
    flush();
    // Collapse and open the structure: its rows come from the incremental flattening.
    t.structure.toggleExpanded();
    t.change('s');
    flush();
    t.structure.toggleExpanded();
    t.change('s');
    flush();
    const host = t.container.querySelector('[data-scene-row="x.row.membrane"]');
    expect(host).not.toBeNull();
    // A full flattening: the row made by the reflow is the same row.
    t.change('s');
    t.change('s');
    t.skip();
    flush();
    const after = t.container.querySelector('[data-scene-row="x.row.membrane"]');
    expect(after).toBe(host);
    t.dispose();
  });
});

describe('the scene tree while searching', () => {
  it('searches again when the scene changes, once for a burst', () => {
    const t = mount({ count: 3, filter: 'lig' });
    expect(t.labels()).toEqual([]);
    const ligand = node('lig', 'component', 'Ligand');
    t.structure.setChildren([...t.components, ligand]);
    t.change('s');
    expect(t.labels()).toEqual(['1ABC', 'Ligand']);
    // A rename that matches: the label is the node's own state, the same node (a new child array is the change).
    const other = node('o', 'component', 'Other');
    t.structure.setChildren([ligand, other]);
    t.change('s');
    expect(t.labels()).toEqual(['1ABC', 'Ligand']);
    t.dispose();
  });
});

describe('an extension row whose view renders later', () => {
  it('gets its level when the tree item appears', async () => {
    const t = mount({ count: 2 });
    const ready = observableState(false);
    t.port.register('x', {
      rows: [
        {
          id: 'x.row.late',
          when: (n) => n.kind === 'structure',
          render: domView(() => (
            <Show when={ready.asReadonly}>{() => <SceneRow label="Late" body={() => <p>Planes</p>} />}</Show>
          )),
        },
      ],
    });
    flush();
    const item = () => t.container.querySelector<HTMLElement>('[data-scene-row="x.row.late"] [role=treeitem]');
    expect(item()).toBeNull();
    ready.set(true);
    flush();
    // The observer of the row's host runs after the render.
    await Promise.resolve();
    expect(item()).not.toBeNull();
    expect(item()!.getAttribute('aria-level')).toBe('2');
    expect(item()!.hasAttribute('aria-expanded')).toBe(false);
    t.dispose();
  });
});

describe('the search field of the Structure tab', () => {
  /** The timers left after the tab is paused, with the search typed or not. */
  const pausedTimers = (typed: boolean) => {
    const scene = mockScene();
    const container = document.createElement('div');
    document.body.append(container);
    const parent = new Scope();
    const dispose = render(() => <StructureTab shell={new ShellModel()} scene={scene} />, container, parent);
    container.querySelector<HTMLButtonElement>('button[data-tip="Search the scene"]')!.click();
    flush();
    const input = container.querySelector<HTMLInputElement>('input')!;
    const all = shownInOrder(container, '[role=treeitem]').length;
    if (typed) {
      input.value = 'zzzz-no-match';
      input.dispatchEvent(new Event('input', { bubbles: true }));
      flush();
    }
    parent.setActive(false);
    const left = vi.getTimerCount();
    vi.advanceTimersByTime(500);
    parent.setActive(true);
    flush();
    vi.advanceTimersByTime(500);
    flush();
    const rows = shownInOrder(container, '[role=treeitem]').length;
    dispose();
    return { left, rows, all };
  };

  it('runs no debounce timer while the tab is paused, and searches what was typed when it shows again', () => {
    vi.useFakeTimers();
    const idle = pausedTimers(false);
    const typed = pausedTimers(true);
    expect(typed.left).toBe(idle.left);
    // Searched once shown again: nothing matches.
    expect(typed.rows).toBeLessThan(typed.all);
  });
});

describe("an extension row's whenChanged", () => {
  it('asks `when` again without a change of the scene: the row comes and goes under its node', () => {
    const t = mount({ count: 3 });
    // Whether the structure has the extension's data (a PAE fetched): a state cell, not a change of the scene.
    const has = observableState(false);
    const told = observableState(0);
    t.port.register('x', {
      rows: [
        {
          id: 'x.row.pae',
          when: (n) => n.kind === 'structure' && has.value,
          whenChanged: told.asReadonly,
          render: domView(() => <SceneRow label="Predicted aligned error" />),
        },
      ],
    });
    const host = () => t.container.querySelector('[data-scene-row="x.row.pae"]');
    expect(t.labels()).not.toContain('Predicted aligned error');
    expect(host()).toBeNull();
    has.set(true);
    told.set(1);
    expect(t.labels()).toContain('Predicted aligned error');
    has.set(false);
    told.set(2);
    expect(t.labels()).not.toContain('Predicted aligned error');
    t.dispose();
  });
});
