/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

/** One thing to do with the current selection (focus it, find it in the tree, make a component of it, …). */
export interface SelectionAction {
  id: string;
  label: string;
  /** One line under the label in the menu: what it does. */
  hint?: string;
  icon: IconName;
  /** Its place among the menu's rows (the viewer's are 0–3, More in the Selection panel last). */
  order?: number;
  run(): void;
}

/**
 * The current selection, for the selection bar (implemented over the structure selection manager): what is
 * selected, and what can be done with it, the menu the summary opens: the viewer's actions and the extensions'
 * (`ui.selectionActions`), by order.
 */
export interface SelectionModel {
  /** "12 residues"; undefined when nothing is selected. */
  summary: ReadonlyState<string | undefined>;
  actions: ReadonlyState<readonly SelectionAction[]>;
  clear(): void;
  /** Something is focused: the Unfocus button beside the bar (or beside the canvas toolbar with nothing selected) */
  focused?: ReadonlyState<boolean>;
  /** Clears the focus, the selection kept (`focus.clear`) */
  unfocus?(): void;
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
