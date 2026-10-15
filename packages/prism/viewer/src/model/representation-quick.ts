/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

// A structure representation's quick settings (parity.md §2): what people change most, above its whole form: the type,
// the color theme (one color picked from a row), the opacity and the size. Drags are one undo step, as a map's levels.

export interface QuickChoice {
  readonly id: string;
  readonly label: string;
}

export interface RepresentationQuickChange {
  readonly type?: string;
  readonly theme?: string;
  /** The one color (the uniform theme). */
  readonly color?: number;
  readonly opacity?: number;
  /** A multiple of the type's default size. */
  readonly size?: number;
}

export interface RepresentationQuickModel {
  /** The common types that apply to its structure (the current one among them). */
  readonly types: ReadonlyState<readonly QuickChoice[]>;
  /** The common color themes that apply (the current one among them). */
  readonly themes: ReadonlyState<readonly QuickChoice[]>;
  readonly type: ReadonlyState<string>;
  readonly theme: ReadonlyState<string>;
  /** Its color when it is one (the uniform theme). */
  readonly color: ReadonlyState<number | undefined>;
  /** 0 to 1; undefined: the type has none. */
  readonly opacity: ReadonlyState<number | undefined>;
  /** A multiple of the type's default size; undefined: the type has none. */
  readonly size: ReadonlyState<number | undefined>;
  /** Changes in quick succession (a slider dragged) are one undo step. */
  change(change: RepresentationQuickChange): void;
  settle(): void;
}
