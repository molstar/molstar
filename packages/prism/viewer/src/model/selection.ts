/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

/** One thing to do with the current selection (color it, make a component of it, focus it, …). */
export interface SelectionAction {
  id: string;
  label: string;
  icon: IconName;
  run(): void;
}

/**
 * The current selection, for the selection bar (implemented over the structure selection manager): what is
 * selected, and what can be done with it. The actions are commands of the `selection` group, so the command menu offers
 * the same ones.
 */
export interface SelectionModel {
  /** "12 residues"; undefined when nothing is selected. */
  summary: ReadonlyState<string | undefined>;
  actions: readonly SelectionAction[];
  clear(): void;
}

/** One structure the selection is in: its name, and where in it ("chain A · 60 atoms"). */
export interface SelectionPart {
  id: string;
  label: string;
  detail: string;
}

/**
 * The selection at the top of the Selection tab (what its sections act on): the summary in Mol*'s words, a line per
 * structure it is in, and Focus and Clear. Over the interaction source and the API (`focus.set`, `select.clear`).
 */
export interface SelectionSummaryModel {
  /** "8 residues"; undefined when nothing is selected. */
  summary: ReadonlyState<string | undefined>;
  parts: ReadonlyState<readonly SelectionPart[]>;
  /** Clicks on the canvas select (the canvas's Select mode); with nothing selected, else the summary offers it. */
  selectMode: ReadonlyState<boolean>;
  /** What a click picks, for people ("Residue", "Chain"). */
  pickLevel: ReadonlyState<string>;
  /** Turns the canvas's Select mode on (`select.mode`). */
  startSelecting(): void;
  focus(): void;
  clear(): void;
}
