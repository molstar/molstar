/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

// The look (api.md §2, "Look"): how representations are drawn. The scene defaults (the Structure panel's sliders) and
// a node's own values (its inspector's sliders), and the overrides between them, marked where they are.

/** One look value some representations keep as their own. */
export interface LookOverrideSummary {
  /** The value's name (`look.set`'s: `material`, `quality`). */
  readonly key: string;
  /** For people: "Material". */
  readonly label: string;
  /** The representations with their own value, by name ("Cartoon of Polymer (1CBS)"). */
  readonly names: readonly string[];
}

/** The scene defaults: the Structure panel's sliders popover. */
export interface SceneLookModel {
  /** The defaults, each change a `look.set` of the scene. */
  readonly form: ParamsFormModel<any>;
  /** The values some representations keep as their own (empty: every one follows the defaults). */
  readonly overrides: ReadonlyState<readonly LookOverrideSummary[]>;
  /** Writes the default of `key` to every representation, overrides included (`look.reset`). */
  applyToAll(key: string): void;
}

/** A structure's, a component's or a representation's look: its inspector's sliders popover. */
export interface NodeLookModel {
  /** The values of the representations under the node (the first one's where they differ); a change sets them all. */
  readonly form: ParamsFormModel<any>;
  /** The values that are not the scene defaults here (on some of the representations under the node). */
  readonly own: ReadonlyState<readonly LookOverrideSummary[]>;
  /** Back to the scene defaults: one value, or all of them. */
  reset(key?: string): void;
  /** Reads the representations again (the scene changed). */
  refresh(): void;
}

/** A tree row's mark: the row has its own look (`own`), or rows under it do (`inside`). */
export interface SceneLookMark {
  readonly kind: 'own' | 'inside';
  /** The tooltip: what is its own ("Its own material and quality"). */
  readonly text: string;
}
