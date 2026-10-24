/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui/dom/view';
import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush } from '@molstar/prism-ui/dom/scope';
import type { SceneNode } from '../../../model/scene.js';
import { mockScene } from '../../../testing/scene.js';
import { extensionUIPort } from '../../../ui/extensions.js';
import type { ViewerUIContext } from '../../../spec.js';
import { SceneRow, SceneRowPickEvent } from '../scene-row.js';
import { SceneTree } from '../scene-tree.js';

// Extensions' rows in the scene tree (viewer-spec.md §3, `rows`): registered through the extension port, shown under
// the nodes their `when` holds for, drawn by `SceneRow` (an eye, a description, controls when opened).

afterEach(() => {
  document.body.innerHTML = '';
});

function mount() {
  const scene = mockScene({ rows: [] });
  const shell = new ShellModel();
  const port = extensionUIPort(shell, {} as never, {} as ViewerUIContext);
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <SceneTree scene={scene} shell={shell} />, container);
  const rowOf = (label: string) =>
    [...container.querySelectorAll('[data-tree-row]')].find(
      (r) => r.querySelector('span[id]')?.textContent === label,
    ) as HTMLElement | undefined;
  return { scene, shell, port, container, dispose, rowOf };
}

describe('extension rows in the scene tree', () => {
  it('shows a row under each structure, with the node as its subject, and removes it with its extension', () => {
    const { scene, port, container, rowOf } = mount();
    const subjects: string[] = [];
    const unregister = port.register('x', {
      rows: [
        {
          id: 'x.row.thing',
          when: (node) => node.kind === 'structure',
          render: domView(({ ctx }) => {
            const node = ctx.subject as SceneNode;
            subjects.push(node.id);
            return <SceneRow label="Thing" description={observableState(`of ${node.id}`).asReadonly} />;
          }),
        },
      ],
    });
    flush();
    const structures = scene.roots.value.filter((n) => n.kind === 'structure');
    expect(structures.length).toBeGreaterThan(0);
    // Under each open structure; not under the other kinds.
    expect(container.querySelectorAll('[data-scene-row="x.row.thing"]').length).toBe(
      structures.filter((s) => s.isExpanded.value).length,
    );
    expect(subjects).toEqual(structures.filter((s) => s.isExpanded.value).map((s) => s.id));
    expect(rowOf('Thing')?.textContent).toContain(`of ${subjects[0]}`);
    // A row is a tree item at the depth of the node's children.
    const host = container.querySelector('[data-scene-row="x.row.thing"]') as HTMLElement;
    expect(host.style.getPropertyValue('--msp-tree-depth')).toBe('1');
    expect(host.querySelector('[role=treeitem]')).not.toBeNull();

    unregister();
    flush();
    expect(container.querySelector('[data-scene-row]')).toBeNull();
  });

  it('opens a closed structure to show it: a node with only rows can be expanded', () => {
    const { scene, port, container } = mount();
    const s = scene.roots.value.find((n) => n.kind === 'structure')!;
    if (s.isExpanded.value) s.toggleExpanded();
    port.register('x', {
      rows: [{ id: 'x.row.a', when: (n) => n === s, render: domView(() => <SceneRow label="A" />) }],
    });
    flush();
    expect(container.querySelector(`[data-scene-row="x.row.a"]`)).toBeNull();
    s.toggleExpanded();
    flush();
    expect(container.querySelectorAll('[data-scene-row="x.row.a"]').length).toBe(1);
  });
});

describe('SceneRow', () => {
  it('toggles through its eye, and shows its controls when opened', () => {
    const shown = observableState(false);
    const onToggle = vi.fn(() => shown.set(!shown.value));
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => (
        <SceneRow
          label="Membrane"
          description={observableState('Off').asReadonly}
          shown={shown.asReadonly}
          onToggle={onToggle}
          body={() => <p class="controls">Planes</p>}
        />
      ),
      container,
    );
    const item = container.querySelector('[role=treeitem]') as HTMLElement;
    const eye = container.querySelector('button[data-tip]') as HTMLButtonElement;
    expect(eye.dataset.tip).toBe('Show');
    expect(item.querySelector('[data-tree-row]')!.className).toMatch(/hidden/);
    eye.click();
    expect(onToggle).toHaveBeenCalledOnce();
    expect(eye.dataset.tip).toBe('Hide');

    expect(container.querySelector('.controls')).toBeNull();
    expect(item.getAttribute('aria-expanded')).toBe('false');
    (container.querySelector('[data-tree-row]') as HTMLElement).click();
    expect(item.getAttribute('aria-expanded')).toBe('true');
    expect(container.querySelector('.controls')?.textContent).toBe('Planes');
  });

  it('is a tree of its own with rows under it: its chevron and a click fold them, its depth indents it', () => {
    const container = document.createElement('div');
    document.body.append(container);
    render(
      () => <SceneRow label="Kinemage" icon="shapes" depth={0} children={() => <SceneRow label="Group" depth={1} />} />,
      container,
    );
    const items = () => [...container.querySelectorAll('[role=treeitem]')] as HTMLElement[];
    const [root] = items();
    expect(root.getAttribute('aria-expanded')).toBe('false');
    expect(items()).toHaveLength(1);
    const expander = root.querySelector('[data-tree-toggle]') as HTMLElement;
    expect(expander.className).not.toMatch(/leaf/);
    (root.querySelector('[data-tree-row]') as HTMLElement).click();
    flush();
    expect(root.getAttribute('aria-expanded')).toBe('true');
    expect(items().map((i) => i.querySelector('span[id]')?.textContent)).toEqual(['Kinemage', 'Group']);
    expect(items()[1].style.getPropertyValue('--msp-tree-depth')).toBe('1');
    // The child is a leaf: no chevron to press.
    expect(items()[1].querySelector('[data-tree-toggle]')!.className).toMatch(/leaf/);
    expander.click();
    flush();
    expect(root.getAttribute('aria-expanded')).toBe('false');
  });

  it('acts on a click with onActivate, in place of being picked', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const onActivate = vi.fn();
    const picked = vi.fn();
    container.addEventListener(SceneRowPickEvent, picked);
    render(() => <SceneRow label="Overview" depth={2} onActivate={onActivate} />, container);
    (container.querySelector('[data-tree-row]') as HTMLElement).click();
    expect(onActivate).toHaveBeenCalledOnce();
    expect(picked).not.toHaveBeenCalled();
  });
});

describe('eyes in the scene tree', () => {
  for (const recycle of [true, false]) {
    it(`draws no eye on a row that shows or hides nothing (Entities, Other), its place kept (${recycle ? 'recycled' : 'keyed'} rows)`, () => {
      const scene = mockScene({ rows: [] });
      const container = document.createElement('div');
      document.body.append(container);
      const dispose = render(() => <SceneTree scene={scene} recycle={recycle} />, container);
      flush();
      const eyeOf = (label: string) => {
        const row = [...container.querySelectorAll('[data-tree-row]')].find(
          (r) => r.querySelector('span[id]')?.textContent === label,
        ) as HTMLElement;
        return row.querySelector('button[aria-label=Hide], button[aria-label=Show]') as HTMLButtonElement;
      };
      for (const label of ['Entities', 'Other']) {
        expect(eyeOf(label).style.visibility).toBe('hidden');
        expect(eyeOf(label).getAttribute('tabindex')).toBe('-1');
      }
      for (const label of ['Components', 'Polymer']) {
        expect(eyeOf(label).style.visibility).toBe('');
        expect(eyeOf(label).hasAttribute('tabindex')).toBe(false);
      }
      dispose();
    });
  }
});
