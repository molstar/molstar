/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

/** One structure's part in a superposition: what of it is used ("1W0E · Chain A"). */
export interface SuperpositionEntry {
  id: string;
  label: string;
}

/**
 * Superposing structures (Mol*'s superposition controls, implemented over the selection and the superposition utils):
 * a method (chains, TM-align, UniProt, ligands, atoms) with its options, the selections it would use (one per structure;
 * UniProt uses every structure), and superpose. Only with two or more structures.
 */
export interface SuperpositionModel {
  /** Two or more structures are loaded. */
  available: ReadonlyState<boolean>;
  /** The method and its options (sequence alignment, trace only). */
  form: ParamsFormModel<any>;
  /** What the current method would superpose. */
  entries: ReadonlyState<readonly SuperpositionEntry[]>;
  /** Enough to superpose; otherwise `hint` says what to select. */
  ready: ReadonlyState<boolean>;
  hint: ReadonlyState<string>;
  focus(id: string): void;
  superpose(): void;
}
