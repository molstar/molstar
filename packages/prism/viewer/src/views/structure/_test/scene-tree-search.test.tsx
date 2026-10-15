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
import { afterEach, describe, expect, it } from 'vitest';
import type { SceneModel, SceneNode, SceneNodeKind } from '../../../model/scene.js';
import { SceneTree } from '../scene-tree.js';

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
  const labels = () =>
    [...container.querySelectorAll('[data-tree-row]')].map((r) => r.querySelector('[id]')!.textContent);
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
    t.contents.toggleExpanded();
    t.contents.toggleExpanded();
    expect(t.labels().filter((l) => l?.startsWith('RIBOSOMAL'))).toHaveLength(5);
  });

  it('searches every row, opening the rows the matches are in', () => {
    const t = mount();
    t.contents.toggleExpanded();
    t.filter.set('l33');
    expect(t.labels()).toEqual(['1JJ2', 'Contents', 'RIBOSOMAL PROTEIN L33']);
    t.filter.set('nothing like it');
    expect(t.container.textContent).toContain('Nothing matches.');
    t.filter.set('');
    expect(t.labels()).toEqual(['1JJ2', 'Contents', 'Polymer']);
  });
});
