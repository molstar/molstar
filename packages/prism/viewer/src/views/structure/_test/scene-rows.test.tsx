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
import type { SceneNode } from '../../../model/scene.js';
import { mockScene } from '../../../testing/scene.js';
import { extensionUIPort } from '../../../ui/extensions.js';
import type { ViewerUIContext } from '../../../spec.js';
import { SceneRow } from '../scene-row.js';
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
    expect(container.querySelector('[data-scene-row]')).toBeNull();
  });

  it('opens a closed structure to show it: a node with only rows can be expanded', () => {
    const { scene, port, container } = mount();
    const s = scene.roots.value.find((n) => n.kind === 'structure')!;
    if (s.isExpanded.value) s.toggleExpanded();
    port.register('x', {
      rows: [{ id: 'x.row.a', when: (n) => n === s, render: domView(() => <SceneRow label="A" />) }],
    });
    expect(container.querySelector(`[data-scene-row="x.row.a"]`)).toBeNull();
    s.toggleExpanded();
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
    const eye = container.querySelector('button[title]') as HTMLButtonElement;
    expect(eye.title).toBe('Show');
    expect(item.querySelector('[data-tree-row]')!.className).toMatch(/hidden/);
    eye.click();
    expect(onToggle).toHaveBeenCalledOnce();
    expect(eye.title).toBe('Hide');

    expect(container.querySelector('.controls')).toBeNull();
    expect(item.getAttribute('aria-expanded')).toBe('false');
    (container.querySelector('[data-tree-row]') as HTMLElement).click();
    expect(item.getAttribute('aria-expanded')).toBe('true');
    expect(container.querySelector('.controls')?.textContent).toBe('Planes');
  });
});
