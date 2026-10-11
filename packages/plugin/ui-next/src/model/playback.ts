/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '../reactive/readonly-state.js';

/** The snapshot stepper's view of the plugin's snapshots (implemented over the snapshot manager). */
export interface SnapshotStepperModel {
  count: ReadonlyState<number>;
  /** 0-based; -1 when no snapshot is current. */
  index: ReadonlyState<number>;
  title: ReadonlyState<string | undefined>;
  /** Markdown, with the snapshot links (`#key`) and markdown extensions of the plugin. */
  description: ReadonlyState<string | undefined>;
  isPlaying: ReadonlyState<boolean>;
  previous(): void;
  next(): void;
  togglePlay(): void;
  /** Follows a link in the description: `#key` applies a snapshot; others are commands or URLs. */
  openLink(href: string): void;
}

/** Model-index controls for a multi-model trajectory (implemented over `UpdateTrajectory`). */
export interface TrajectoryControlsModel {
  count: ReadonlyState<number>;
  /** 0-based. */
  index: ReadonlyState<number>;
  isPlaying: ReadonlyState<boolean>;
  first(): void;
  previous(): void;
  next(): void;
  last(): void;
  seek(index: number): void;
  togglePlay(): void;
}

/** The label under the cursor (implemented over the plugin's hover labels). */
export interface HoverLabelModel {
  /** Lines to show (markdown-free text); empty when nothing is hovered. */
  lines: ReadonlyState<readonly string[]>;
}
