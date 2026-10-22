/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '@molstar/prism-ui/components/icon/icon';
import type { CommandValues } from '@molstar/prism-ui/core';

// The command menu's own forms for a component's look (api.md §9, "Own forms"): Show as… picks a representation type
// among those the component supports, with the type's main options; Color by theme… picks a color theme among those
// that apply, drawn with its colors, and the uniform theme's color. Read for the values the form holds now (the
// component and structure the inspector row filled in), so the view knows nothing of the plugin.

export interface StyleChoice {
  readonly id: string;
  readonly label: string;
  /** A type's icon (as the inspector's representation rows draw it). */
  readonly icon?: IconName;
  /** A theme's colors on the structure, as a CSS background (from its legend); none for one color. */
  readonly swatch?: string;
  /** Its group (a theme's registry category: Atom Property, Chain Property, …); listed under it. */
  readonly category?: string;
  /** The component shows this type already (adding it again updates that one). */
  readonly shown?: boolean;
}

/** One of a type's main options, as the form draws it. */
export type TypeOption =
  | {
      readonly kind: 'range';
      readonly name: string;
      readonly label: string;
      readonly min: number;
      readonly max: number;
      readonly step: number;
      readonly default: number;
    }
  | {
      readonly kind: 'choice';
      readonly name: string;
      readonly label: string;
      readonly options: readonly (readonly [string, string])[];
      readonly default: string;
    };

export interface StyleFormsModel {
  /**
   * The representation types the component the values name supports (those it shows marked); with none named, those
   * its structure (the only one) supports; empty when the values name neither.
   */
  types(values: CommandValues): readonly StyleChoice[];
  /** A type's main options (size, quality, opacity) on that component, at their defaults there. */
  typeOptions(type: string, values: CommandValues): readonly TypeOption[];
  /** The color themes that apply to the structure the values name, in the registry's order, by category. */
  themes(values: CommandValues): readonly StyleChoice[];
}
