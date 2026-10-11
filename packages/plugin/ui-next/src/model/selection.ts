/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

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
