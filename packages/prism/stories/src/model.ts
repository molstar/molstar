/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui';
import type { SnapshotStepperModel } from '@molstar/prism-ui/model/playback';

// The stories' view-model contracts (an app of their own, not a mode: the Prism viewer's `_spec/README.md`, decision
// 3): a story is a document (scenes described by schemas and MolViewSpec), the canvas renders it; editing a story edits
// the document, never the plugin's state directly, so a story stays serializable and reproducible. A first sketch for
// the mocks: the scenes and the player; the rest of the story UX (scene editing, assets, the narrative text,
// publishing) is still to be designed.

/** One scene of a story: what the canvas shows (a MolViewSpec snapshot, by key) and what it says. */
export interface StoryScene {
  id: string;
  title: string;
  /** Markdown. */
  caption?: string;
  /** How long it plays (ms); unset: until the reader moves on. */
  duration?: number;
}

export interface StoryModel {
  title: ReadonlyState<string>;
  scenes: ReadonlyState<readonly StoryScene[]>;
  /** The scene being edited and shown. */
  current: ReadonlyState<string | undefined>;
  select(id: string): void;
  /** A new scene after the current one, from what the canvas shows. */
  addScene(): void;
  removeScene(id: string): void;
  /** The scenes as the canvas plays them (the playback card). */
  player: SnapshotStepperModel;
}
