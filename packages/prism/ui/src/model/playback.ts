/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MarkdownHost } from '../components/markdown/markdown.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

/** The snapshot stepper's view of the plugin's snapshots (implemented over the snapshot manager). */
export interface SnapshotStepperModel {
  count: ReadonlyState<number>;
  /** 0-based; -1 when no snapshot is current. */
  index: ReadonlyState<number>;
  title: ReadonlyState<string | undefined>;
  /** Every snapshot's title, in order (the position dots' tooltips, the list when there are many). */
  titles: ReadonlyState<readonly string[]>;
  /** Markdown, with the snapshot links (`#key`) and markdown extensions of the plugin. */
  description: ReadonlyState<string | undefined>;
  /** Whether the description is plain text (a snapshot's `descriptionFormat`), shown as it is. */
  plainText?: ReadonlyState<boolean>;
  /** The page's part of the description: command links, images, the audio player (Mol*'s markdown extensions). */
  markdown?: MarkdownHost;
  isPlaying: ReadonlyState<boolean>;
  previous(): void;
  next(): void;
  goTo(index: number): void;
  togglePlay(): void;
  /** Follows a link in the description: `#key` applies a snapshot; others are commands or URLs. */
  openLink(href: string): void;
}

/**
 * The current snapshot's transition (frames played after it is applied; v5's snapshot animation slider): where it is,
 * a scrubber through it, play from the start and stop. Its length is 0 when the snapshot has none.
 */
export interface SnapshotTransitionModel {
  /** ms; 0: no transition. */
  duration: ReadonlyState<number>;
  /** The shortest frame (ms): the scrubber's step. */
  step: ReadonlyState<number>;
  /** Where it is (ms). */
  time: ReadonlyState<number>;
  isPlaying: ReadonlyState<boolean>;
  seek(time: number): void;
  /** Plays it from its start, or stops it. */
  togglePlay(): void;
}

/**
 * Model-index controls for a multi-model trajectory (implemented over `UpdateTrajectory`). A trajectory belongs to a
 * structure (its frame is the structure's model index): with several loaded, the bar follows the focused one.
 */
export interface TrajectoryControlsModel {
  /** The structure it belongs to ("1TQN"). */
  label: ReadonlyState<string>;
  count: ReadonlyState<number>;
  /** 0-based. */
  index: ReadonlyState<number>;
  isPlaying: ReadonlyState<boolean>;
  seek(index: number): void;
  togglePlay(): void;
}

/** The label under the cursor (implemented over the plugin's hover labels). */
export interface HoverLabelModel {
  /** Lines to show (markdown-free text); empty when nothing is hovered. */
  lines: ReadonlyState<readonly string[]>;
}
