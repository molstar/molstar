/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ApiFailure } from './area.js';

// The target resolver operations get (targets.md): a target, in any of its forms, to loci on a structure and "what it
// selected". Built over the plugin in plan phase 2.1, on the `targets/` port (track A8).

export interface TargetResolver {
  /** Resolves a target on a structure (a state ref or a label; absent: the structure the person works with). */
  resolve(target: unknown, structure?: string): Promise<ResolvedTarget | ApiFailure>;
}

export interface ResolvedTarget {
  readonly ok: true;
  /** The structure it was resolved on: its state ref and label. */
  readonly structure: { readonly ref: string; readonly label: string };
  /** The loci (a `StructureElement.Loci`), for the plugin's managers. */
  readonly loci: unknown;
  /** What it selected, in words ("42 residues in chains A and B"). */
  readonly summary: string;
  /** The target as stored: MolQL JSON. */
  readonly stored: unknown;
}

/** Until the resolver is built: says so. */
export const NoTargets: TargetResolver = {
  async resolve() {
    return { ok: false, code: 'unavailable', message: 'Targets are not available in this viewer yet.' };
  },
};
