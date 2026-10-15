/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

/** What a style layer row shows: the atoms it is on, and its value. */
export interface StyleLayerView {
  /** A color layer's color (`#rrggbb`): the swatch; else the kind's icon. */
  readonly swatch?: string;
  readonly icon: IconName;
  /** Where it is ("A:10-17"; "1CBS · A:10-17" when several structures have layers). */
  readonly where: string;
  /** Its value, short ("50% transparent", "glow 0.5"); empty for a color (the swatch says it). */
  readonly value: string;
  /** All of it in words, for the tooltip ("Color #e0523a on 8 residues (60 atoms) in chain A of 1CBS"). */
  readonly detail: string;
}

/** One applied layer (Mol*'s overpaint, transparency, emissive, substance or clipping on part of a structure). */
export interface StyleLayerItem {
  /** Stable while the row lives: an edit of its value keeps the row (and its open editor). */
  readonly key: string;
  readonly view: ReadonlyState<StyleLayerView>;
  /** Its value can be changed in place (a color, a transparency, a glow; not a reset, a material or clipping). */
  readonly editable: boolean;
  /** Its value in a form (a color, an amount), loaded with the current value; undefined: not changed in place. */
  editor(): ParamsFormModel<any> | undefined;
  remove(): void;
}

/**
 * The style layers applied (the Selection tab's "Color, transparency, glow" section, under its form): a row each, with
 * an inline edit of its value and remove, and Clear all. The rows follow the state while they show.
 */
export interface StyleLayersModel {
  readonly items: ReadonlyState<readonly StyleLayerItem[]>;
  clear(): void;
}
