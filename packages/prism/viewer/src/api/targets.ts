/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureElement } from '@molstar/model/model/structure';
import type { QueryLanguage, Selected } from '../targets/types.js';
import type { ApiFailure } from './area.js';

// The target resolver operations get (targets.md): a target, in any of its forms, to loci on the structures it applies
// to, with "what it selected". Built over the plugin in `sources/targets.ts`, on the `targets/` port.

export interface ResolveOptions {
  /** A structure by state ref or label; absent: the one the person works with (targets.md §1). */
  readonly structure?: string;
  /** How a query text reads; PyMOL by default. */
  readonly language?: QueryLanguage;
}

/** The target on one structure. */
export interface ResolvedPart {
  readonly structure: { readonly ref: string; readonly label: string };
  readonly loci: StructureElement.Loci;
  readonly selected: Selected;
}

export interface ResolvedTarget {
  readonly ok: true;
  /** One per structure it was resolved on. */
  readonly parts: readonly ResolvedPart[];
  /** Which rule chose the structures: named, the focused one, the selected, the only one, or every loaded one. */
  readonly scope: 'named' | 'focused' | 'selected' | 'only' | 'all';
  /** What it selected, in words ("42 residues (350 atoms) in chains A and B"). */
  readonly summary: string;
  /** It selects nothing anywhere. */
  readonly empty: boolean;
  /** The target as stored: MolQL JSON (`{ molql }`). */
  readonly stored: unknown;
}

export interface TargetResolver {
  resolve(target: unknown, options?: ResolveOptions): Promise<ResolvedTarget | ApiFailure>;
}

/** Without a plugin's structures (tests of the frame): says so. */
export const NoTargets: TargetResolver = {
  async resolve() {
    return { ok: false, code: 'unavailable', message: 'Targets are not available in this viewer.' };
  },
};
