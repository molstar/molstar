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
  /** What its frames are called (an ensemble's models, a simulation's frames); `models` when not given. */
  unit?: ReadonlyState<'models' | 'frames'>;
  count: ReadonlyState<number>;
  /** 0-based. */
  index: ReadonlyState<number>;
  isPlaying: ReadonlyState<boolean>;
  seek(index: number): void;
  togglePlay(): void;
  /** A simulation's periodic box: how its frames show it, when the structure is in one (the bar's box menu). */
  periodic?: PeriodicBoxModel;
}

/** How a simulation box's frames are shown: molecules whole, in the box (centered on the polymers), without jumps. */
export type PeriodicBoxMode = 'whole' | 'box' | 'box-centered' | 'nojump' | 'off';

export interface PeriodicBoxModel {
  /** The mode now; absent: the structure is not in a simulation box (no menu). */
  mode: ReadonlyState<PeriodicBoxMode | undefined>;
  set(mode: PeriodicBoxMode): void;
  /** The neighbouring boxes (periodic images): shown or not, and the switch. */
  images?: { shown: ReadonlyState<boolean>; set(shown: boolean): void };
}

/** A run of a hover label's line: text with its emphasis (never markup: the view builds elements from these). */
export interface HoverRun {
  readonly text: string;
  readonly bold?: boolean;
  readonly italic?: boolean;
  readonly small?: boolean;
  readonly code?: boolean;
  readonly sub?: boolean;
  readonly sup?: boolean;
  /** Link text (a label takes no pointer: shown as a link, not followed). */
  readonly link?: boolean;
}

/** A formatted line: its runs and its role (a title, text, a list item). */
export interface HoverRichLine {
  readonly runs: readonly HoverRun[];
  readonly role: 'title' | 'text' | 'item';
}

/** A hover label's line: plain text (the first is the title), or formatted runs. */
export type HoverLine = string | HoverRichLine;

/** The label under the cursor (implemented over the plugin's hover labels). */
export interface HoverLabelModel {
  /** Lines to show; empty when nothing is hovered. */
  lines: ReadonlyState<readonly HoverLine[]>;
}
