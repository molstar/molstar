/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { StructureElement } from '@molstar/model/model/structure';
import type { Target } from '../targets/types.js';

// The `target` param (targets.md §1): a keyword, an address list, a query text or a selector (schema items, `{ molql }`),
// and, from code in the page (the sequence link), loci directly (`{ loci }`, never stored). A plain value param: its
// command form shows a text field (`commands/operations.ts`), its tool field is the same string (agent.md §3).

/** What an operation's `target` takes. */
export type TargetValue = Target | { readonly loci: StructureElement.Loci };

export function TargetParam(description = 'What to act on: ligand, A:120-130, resn HEM, …'): PD.Value<TargetValue> {
  return PD.Value<TargetValue>('', { description });
}

/** The query language of a target's text (targets.md §1): PyMOL unless said. */
export function LanguageParam() {
  return PD.Select<'pymol' | 'vmd' | 'jmol' | 'molscript'>('pymol', [
    ['pymol', 'PyMOL'],
    ['vmd', 'VMD'],
    ['jmol', 'Jmol'],
    ['molscript', 'MolScript'],
  ]);
}

/** Whether a value names loci directly. */
export function isLociTarget(t: unknown): t is { readonly loci: StructureElement.Loci } {
  return !!t && typeof t === 'object' && 'loci' in t && !!(t as { loci: unknown }).loci;
}
