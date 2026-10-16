/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState, type ObservableState } from '@molstar/prism-ui/reactive/observable-state';
import type { SceneModel, SceneNode, SceneNodeKind } from '../model/scene.js';

// Fixtures of the performance stories: a scene of one structure with many components, built from plain states (no
// viewer), the shape the real scene tree reads.

export interface SceneFixtureNode {
  readonly node: SceneNode;
  readonly label: ObservableState<string>;
  readonly expanded: ObservableState<boolean>;
}

export interface SceneFixture {
  scene: SceneModel;
  structure: SceneFixtureNode;
  components: SceneFixtureNode[];
  /** Representation rows under all components (shown when those are open). */
  representations: number;
}

function node(id: string, kind: SceneNodeKind, label: string, children: SceneFixtureNode[] = [], expanded = false) {
  const text = observableState(label);
  const isExpanded = observableState(expanded);
  const hidden = observableState(false);
  const n = {
    id,
    kind,
    icon: kind === 'representation' ? 'cartoon' : 'structure',
    label: text.asReadonly,
    description: observableState<string | undefined>(kind === 'representation' ? 'Cartoon · Chain ID' : '120 atoms')
      .asReadonly,
    isHidden: hidden.asReadonly,
    isExpanded: isExpanded.asReadonly,
    children: observableState<readonly SceneNode[]>(children.map((c) => c.node)).asReadonly,
    toggleVisibility: () => hidden.set(!hidden.value),
    toggleExpanded: () => isExpanded.set(!isExpanded.value),
    highlight: () => {},
    focus: () => {},
    remove: () => {},
  } as unknown as SceneNode;
  return { node: n, label: text, expanded: isExpanded } satisfies SceneFixtureNode;
}

/** One structure with `components` components (closed), each with `representations` representations. */
export function makeScene(components = 1000, representations = 3): SceneFixture {
  const comps = Array.from({ length: components }, (_, i) =>
    node(
      `c${i}`,
      'component',
      `Component ${i}`,
      Array.from({ length: representations }, (_, r) => node(`c${i}r${r}`, 'representation', `Representation ${r}`)),
    ),
  );
  const entities = Array.from({ length: 12 }, (_, i) => node(`e${i}`, 'entity', `Entity ${i}`));
  const contents = node('contents', 'contents', 'Contents', entities);
  const structure = node('s', 'structure', '1ABC', [contents, ...comps], true);
  const picked = observableState<SceneNode | undefined>(undefined);
  const scene = {
    roots: observableState([structure.node]).asReadonly,
    picked: picked.asReadonly,
    pick: (n: SceneNode | undefined) => picked.set(n),
  } as unknown as SceneModel;
  return { scene, structure, components: comps, representations: components * representations };
}
