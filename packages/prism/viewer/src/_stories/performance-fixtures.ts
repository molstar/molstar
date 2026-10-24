/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState, type ObservableState } from '@molstar/prism-ui/reactive/observable-state';
import type { SceneModel, SceneNode, SceneNodeKind } from '../model/scene.js';
import { SceneChangeLog } from '../model/scene-changes.js';

// Fixtures of the performance stories: a scene of one structure with many components, built from plain states (no
// viewer), the shape the real scene tree reads, its structure's changes told as the scene source tells them
// (`SceneModel.changes`: a node opened or closed, by id).

export interface SceneFixtureNode {
  readonly node: SceneNode;
  readonly label: ObservableState<string>;
  readonly expanded: ObservableState<boolean>;
  readonly hidden: ObservableState<boolean>;
  readonly children: ObservableState<readonly SceneNode[]>;
  /** With `swatches` only. */
  readonly swatch?: ObservableState<number | undefined>;
}

export interface SceneFixture {
  scene: SceneModel;
  structure: SceneFixtureNode;
  components: SceneFixtureNode[];
  /** Representation rows under all components (shown when those are open). */
  representations: number;
}

function node(
  id: string,
  kind: SceneNodeKind,
  label: string,
  children: SceneFixtureNode[] = [],
  expanded = false,
  withSwatch = false,
  log?: SceneChangeLog,
) {
  const text = observableState(label);
  const swatch = withSwatch ? observableState<number | undefined>(0x2e86c1) : undefined;
  const isExpanded = observableState(expanded);
  const kids = observableState<readonly SceneNode[]>(children.map((c) => c.node));
  // Opened or closed (by the row, or a test setting it), its children changed: the scene's change stream says so.
  if (log) {
    isExpanded.asReadonly.listen?.(() => log.touch(id));
    kids.asReadonly.listen?.(() => log.touch(id));
  }
  const hidden = observableState(false);
  const n = {
    id,
    kind,
    icon: kind === 'representation' ? 'cartoon' : 'structure',
    label: text.asReadonly,
    description: observableState<string | undefined>(kind === 'representation' ? 'Cartoon · Chain ID' : '120 atoms')
      .asReadonly,
    isHidden: hidden.asReadonly,
    ...(swatch ? { swatch: swatch.asReadonly } : {}),
    isExpanded: isExpanded.asReadonly,
    children: kids.asReadonly,
    toggleVisibility: () => hidden.set(!hidden.value),
    toggleExpanded: () => isExpanded.set(!isExpanded.value),
    highlight: () => {},
    focus: () => {},
    remove: () => {},
  } as unknown as SceneNode;
  return { node: n, label: text, expanded: isExpanded, hidden, children: kids, swatch } satisfies SceneFixtureNode;
}

/**
 * One structure with `components` components (closed), each with `representations` representations; with `swatches`,
 * each component has a color swatch.
 */
export function makeScene(components = 1000, representations = 3, swatches = false): SceneFixture {
  const log = new SceneChangeLog();
  const comps = Array.from({ length: components }, (_, i) =>
    node(
      `c${i}`,
      'component',
      `Component ${i}`,
      Array.from({ length: representations }, (_, r) => node(`c${i}r${r}`, 'representation', `Representation ${r}`)),
      false,
      swatches,
      log,
    ),
  );
  const entities = Array.from({ length: 12 }, (_, i) => node(`e${i}`, 'entity', `Entity ${i}`));
  const contents = node('contents', 'contents', 'Entities', entities, false, false, log);
  const structure = node('s', 'structure', '1ABC', [contents, ...comps], true, false, log);
  const picked = observableState<SceneNode | undefined>(undefined);
  const scene = {
    roots: observableState([structure.node]).asReadonly,
    picked: picked.asReadonly,
    pick: (n: SceneNode | undefined) => picked.set(n),
    changes: log.changes,
  } as unknown as SceneModel;
  return { scene, structure, components: comps, representations: components * representations };
}
