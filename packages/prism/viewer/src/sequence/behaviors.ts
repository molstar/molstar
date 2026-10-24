/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import type { IconName, LegendItem, ReadonlyState } from '@molstar/prism-ui';
import type { Op } from '@molstar/seq-core/format';

// Sequence behaviors (sequence.md §3): small rules that keep their part of the sequence's document in step with what it
// shows. The viewer builds the document's base for the structure it follows (the ruler, the letters, the secondary
// structure); a behavior adds to it by patch once that document is in, with refs that start with its id. A new document
// (another structure, another letter theme) has none of a behavior's nodes: `shown` changes and the behavior adds them
// again. Extensions bring theirs (`ViewerExtension.sequence`): model-archive's pLDDT track (README decision 27).

/** What the sequence shows, once its document is in: a new value for each document. */
export interface SequenceShown {
  /** Tells one structure from another. */
  readonly key: string;
  readonly title: string;
  /** The structure the 3D view draws (its model carries the custom properties). */
  readonly structure: Structure;
  /** The complex's ref: a computed data node names it in `options.of`. */
  readonly complex: string;
  /** Where data nodes go (the document's `data` node). */
  readonly data: string;
  /** The section the tracks are in, and its last track: a behavior's tracks go after it (none: the section's end). */
  readonly section: string;
  readonly after?: string;
}

/** What a track's menu offers: a coloring of the 3D view's polymer, by a color theme's name. */
export interface SequenceTrackMenuItem {
  readonly id: string;
  readonly label: string;
  readonly icon: IconName;
  /** The 3D color theme (`plddt-confidence`); `default` is the viewer's default coloring. */
  readonly theme: string;
}

/** What a behavior adds beside its nodes: its tracks' menus (by track ref) and legend blocks. */
export interface SequenceOffer {
  readonly menus?: Readonly<Record<string, readonly SequenceTrackMenuItem[]>>;
  /** Its tracks' names for people (the Visible tracks menu, sequence.md §8), by ref; absent: the behavior's title. */
  readonly labels?: Readonly<Record<string, string>>;
  readonly legend?: readonly { readonly title: string; readonly items: readonly LegendItem[] }[];
}

export interface SequenceBehaviorContext {
  /** What the sequence shows; undefined: nothing. */
  readonly shown: ReadonlyState<SequenceShown | undefined>;
  /**
   * Patches the document's base (Seq*'s `applyPatch`, all or nothing): for the document `shown` names, so a patch made
   * for one that was replaced meanwhile is dropped (false), as is a refused one.
   */
  patch(shown: SequenceShown, ops: readonly Op[]): Promise<boolean>;
  /** Its menus and legend blocks now (replaces its last offer; undefined withdraws it). */
  offer(offer: SequenceOffer | undefined): void;
}

export interface SequenceBehavior {
  /** `<owner>.sequence.<name>` (viewer-spec.md §4). */
  readonly id: string;
  readonly title: string;
  /** Starts it; returns its dispose (which withdraws its offer; its nodes go with the next document). */
  attach(ctx: SequenceBehaviorContext): () => void;
}
