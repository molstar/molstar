/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { domView } from '@molstar/prism-ui/dom/view';
import { afterEach, describe, expect, it } from 'vitest';
import { flush } from '@molstar/prism-ui/dom/scope';
import type { SceneModel, SceneNode, SceneNodeKind } from '../../../model/scene.js';
import { extensionUIPort } from '../../../ui/extensions.js';
import type { ViewerUIContext } from '../../../spec.js';
import { SceneRow } from '../scene-row.js';
import { SceneTree } from '../scene-tree.js';
import { shownInOrder } from '@molstar/prism-ui/testing/virtual-list';

// A big structure in the tree: Contents shows its first entities and Show all for the rest (the few again once closed
// and opened); a search shows the matching rows with the rows they are in, opened, whatever was shown before.

afterEach(() => {
  document.body.innerHTML = '';
});

function node(id: string, kind: SceneNodeKind, label: string, children: SceneNode[] = [], expanded = false) {
  const isExpanded = observableState(expanded);
  const n: SceneNode = {
    id,
    kind,
    icon: 'structure',
    label: observableState(label).asReadonly,
    description: observableState<string | undefined>(undefined).asReadonly,
    isHidden: observableState(false).asReadonly,
    isExpanded: isExpanded.asReadonly,
    children: observableState<readonly SceneNode[]>(children).asReadonly,
    toggleVisibility: () => {},
    toggleExpanded: () => isExpanded.set(!isExpanded.value),
    highlight: () => {},
    focus: () => {},
    remove: () => {},
  } as unknown as SceneNode;
  return n;
}

function mount() {
  const entities = Array.from({ length: 36 }, (_, i) => node(`e${i}`, 'entity', `RIBOSOMAL PROTEIN L${i}`));
  const contents = node('contents', 'contents', 'Contents', entities, true);
  const structure = node('s', 'structure', '1JJ2', [contents, node('c', 'component', 'Polymer')], true);
  const scene = {
    roots: observableState([structure]).asReadonly,
    picked: observableState<SceneNode | undefined>(undefined).asReadonly,
    pick: () => {},
  } as unknown as SceneModel;
  const filter = observableState('');
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <SceneTree scene={scene} filter={filter.asReadonly} />, container);
  const labels = () => {
    // Read after what was set is rendered (the runtime's flush: the rows are flattened once a burst).
    flush();
    return shownInOrder(container, '[data-tree-row]').map((r) => r.querySelector('[id]')!.textContent);
  };
  return { contents, filter, container, labels };
}

describe('the scene tree with many rows', () => {
  it('shows a few entities, then all on Show all, and the few again once Contents is closed and opened', () => {
    const t = mount();
    expect(t.labels().filter((l) => l?.startsWith('RIBOSOMAL'))).toHaveLength(5);
    const more = t.container.querySelector<HTMLButtonElement>('button[class*="more"]')!;
    expect(more.textContent).toContain('Show all 36');
    more.click();
    expect(t.labels().filter((l) => l?.startsWith('RIBOSOMAL'))).toHaveLength(36);
    // Closed, then opened: two presses (each flushed, as an event's dispatch is).
    t.contents.toggleExpanded();
    flush();
    t.contents.toggleExpanded();
    expect(t.labels().filter((l) => l?.startsWith('RIBOSOMAL'))).toHaveLength(5);
  });

  it('looks into lazy contents once per search, leaving them lazy (later syncs do not compute them)', () => {
    const entities = observableState<readonly SceneNode[]>([]);
    const rows = [node('e1', 'entity', 'RETINOIC ACID BINDING PROTEIN'), node('e2', 'entity', 'WATER')];
    const loads: (boolean | undefined)[] = [];
    const contents = {
      ...node('contents', 'contents', 'Contents'),
      children: entities.asReadonly,
      lazy: true,
      load(peek?: boolean) {
        loads.push(peek);
        if (!peek) contents.lazy = false;
        entities.set(rows);
      },
    } as SceneNode;
    const structure = node('s', 'structure', '1CBS', [contents, node('c', 'component', 'Polymer')], true);
    const scene = {
      roots: observableState([structure]).asReadonly,
      picked: observableState<SceneNode | undefined>(undefined).asReadonly,
      pick: () => {},
    } as unknown as SceneModel;
    const filter = observableState('');
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <SceneTree scene={scene} filter={filter.asReadonly} />, container);
    const labels = () => {
      // Read after what was set is rendered (the runtime's flush: the rows are flattened once a burst).
      flush();
      return shownInOrder(container, '[data-tree-row]').map((r) => r.querySelector('[id]')!.textContent);
    };

    filter.set('retinoic');
    expect(labels()).toEqual(['1CBS', 'Contents', 'RETINOIC ACID BINDING PROTEIN']);
    filter.set('water');
    expect(labels()).toEqual(['1CBS', 'Contents', 'WATER']);
    // Looks only (never a load that keeps them): one walk per search whatever the rows reading it, and the opening's.
    expect(loads).toEqual([true, true, true]);
    expect(contents.lazy).toBe(true);
    filter.set('');
    expect(labels()).toEqual(['1CBS', 'Contents', 'Polymer']);
  });

  it('searches every row, opening the rows the matches are in', () => {
    const t = mount();
    t.contents.toggleExpanded();
    t.filter.set('l33');
    expect(t.labels()).toEqual(['1JJ2', 'Contents', 'RIBOSOMAL PROTEIN L33']);
    t.filter.set('nothing like it');
    flush();
    expect(t.container.textContent).toContain('Nothing matches.');
    t.filter.set('');
    expect(t.labels()).toEqual(['1JJ2', 'Contents', 'Polymer']);
  });
});

describe('the scene tree search and extension rows', () => {
  it('finds a structure by the title of an extension row under it, opened to show the row', () => {
    const structure = node('s', 'structure', '1JJ2', [node('c', 'component', 'Polymer')], false);
    const scene = {
      roots: observableState([structure]).asReadonly,
      picked: observableState<SceneNode | undefined>(undefined).asReadonly,
      pick: () => {},
    } as unknown as SceneModel;
    const shell = new ShellModel();
    const port = extensionUIPort(shell, {} as never, {} as ViewerUIContext);
    port.register('x', {
      rows: [
        {
          id: 'x.row.membrane',
          when: (n) => n.kind === 'structure',
          render: domView(() => <SceneRow label="Membrane" />),
        },
      ],
    });
    const filter = observableState('');
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <SceneTree scene={scene} shell={shell} filter={filter.asReadonly} />, container);
    const labels = () => {
      // Read after what was set is rendered (the runtime's flush: the rows are flattened once a burst).
      flush();
      return shownInOrder(container, '[data-tree-row]').map((r) => r.querySelector('[id]')!.textContent);
    };

    filter.set('membrane');
    flush();
    expect(container.textContent).not.toContain('Nothing matches.');
    expect(labels()).toEqual(['1JJ2', 'Membrane']);
    filter.set('nothing like it');
    flush();
    expect(container.textContent).toContain('Nothing matches.');
  });
});
