/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { CifBlock } from '@molstar/io/reader/cif/data-model';

/** How a tier C component links (plan §5): polymer residues, glycans and covalent ligands are reported apart. */
export type LinkKind = 'polymer' | 'saccharide' | 'ligand';

export interface Representative {
  /** `chem_comp.pdbx_model_coordinates_db_code`, upper case; absent when the CCD names none */
  entry?: string;
  /** `polymer` and `saccharide` from the type alone; `ligand` for the rest, a candidate only when listed as bound */
  kind: LinkKind;
}

export function representative(block: CifBlock): Representative {
  const chemComp = block.categories.chem_comp;
  const code = chemComp?.getField('pdbx_model_coordinates_db_code')?.str(0).trim().toUpperCase();
  const type = chemComp?.getField('type')?.str(0).toUpperCase() ?? '';
  const kind: LinkKind = type.includes('SACCHARIDE') ? 'saccharide' : type.includes('LINKING') ? 'polymer' : 'ligand';
  return { entry: code && code !== '?' && code !== '.' ? code : void 0, kind };
}

/**
 * A tier C candidate (plan §3.2): a polymer or saccharide type, or a component its entry lists as bound
 * (`nonpolymer_bound_components`, which counts metal coordination too: the entry file decides).
 */
export function isCandidate(id: string, rep: Representative, bound: Map<string, string[] | null>) {
  if (!rep.entry) return false;
  if (rep.kind !== 'ligand') return true;
  return bound.get(rep.entry)?.includes(id) ?? false;
}
