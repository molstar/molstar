/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Column } from '@molstar/core/data/db';
import type { CifBlock } from '@molstar/io/reader/cif/data-model';
import { ElementSymbol } from '@molstar/model/model/structure/model/types';
import { AtomNumber } from '@molstar/model/model/structure/model/properties/atomic/measures';
import { isMetal } from '@molstar/model/model/structure/model/properties/atomic/types';
import { isHydrogen } from './input.ts';
import type { Tier } from './tiers.ts';

/** The V2000 counts line has three digits. */
const MaxAtoms = 999;

export interface CcdCheck {
  /** why the component is left out of the tier (plan §5.1), from the CCD alone */
  excluded?: string;
  /** leaving atoms without coordinates (tier B): compared without them (plan §4.2) */
  leavingMissing: number;
}

/**
 * Whether a component is benchmarked in a tier. Every atom written (heavy atoms, and H in tier 0) needs the tier's
 * coordinates, except leaving atoms in model coordinates (OXT, OP3, a linked O1: absent where the residue was linked).
 */
export function checkCcd(block: CifBlock, tier: Tier): CcdCheck {
  const atoms = block.categories.chem_comp_atom;
  if (!atoms || atoms.rowCount === 0) return { excluded: 'no atoms', leavingMissing: 0 };
  const typeSymbol = atoms.getField('type_symbol')!;
  const leaving = atoms.getField('pdbx_leaving_atom_flag');
  // tier C's coordinates are the entry's: checked on the copy (`entry.ts`)
  const x = tier.entry
    ? void 0
    : atoms.getField(tier.coordinates === 'ideal' ? 'pdbx_model_Cartn_x_ideal' : 'model_Cartn_x');
  const exclude = (excluded: string) => ({ excluded, leavingMissing: 0 });

  let heavy = 0;
  let written = 0;
  let leavingMissing = 0;
  for (let i = 0; i < atoms.rowCount; i++) {
    const element = ElementSymbol(typeSymbol.str(i));
    if (isMetal(element)) return exclude('organometallic (metal atom)');
    if (AtomNumber(element) === 0) return exclude(`unknown element ${typeSymbol.str(i)}`);
    const h = isHydrogen(element);
    if (h && !tier.hydrogens) continue;
    if (tier.entry || (x && x.valueKind(i) === Column.ValueKind.Present)) {
      written++;
      if (!h) heavy++;
    } else if (tier.coordinates === 'model' && leaving?.str(i).toUpperCase() === 'Y') {
      leavingMissing++;
    } else {
      return exclude(`${tier.coordinates} coordinates incomplete`);
    }
  }
  if (heavy <= 1) return exclude('single heavy atom');
  const cage = nonClassical(block);
  if (cage) return exclude(`non-classical bonding (${cage})`);
  if (written > MaxAtoms) return exclude(`over ${MaxAtoms} atoms (V2000)`);
  // tier 0 compares with the CCD's own InChI; the perception tiers with RDKit's of the CCD graph (plan §4.2)
  if (!tier.perceive && !ccdInchi(block)) return exclude('no CCD InChI');
  return { leavingMissing };
}

/** Most bonds of an atom in a classical structure: past them, a cage (carborane) or a bridging H (diborane). */
const MaxBonds: Record<string, number> = { H: 1, B: 4, C: 4 };

/** The first atom of the CCD graph with more bonds than a classical structure allows, as `name: count`. */
function nonClassical(block: CifBlock): string | undefined {
  const atoms = block.categories.chem_comp_atom;
  const bonds = block.categories.chem_comp_bond;
  if (!bonds) return void 0;
  const elementOf = new Map<string, string>();
  const atomId = atoms.getField('atom_id')!;
  const typeSymbol = atoms.getField('type_symbol')!;
  for (let i = 0; i < atoms.rowCount; i++) elementOf.set(atomId.str(i), typeSymbol.str(i).toUpperCase());
  const count = new Map<string, number>();
  for (const field of ['atom_id_1', 'atom_id_2']) {
    const f = bonds.getField(field)!;
    for (let i = 0; i < bonds.rowCount; i++) count.set(f.str(i), (count.get(f.str(i)) ?? 0) + 1);
  }
  for (const [name, n] of count) {
    const max = MaxBonds[elementOf.get(name) ?? ''];
    if (max !== void 0 && n > max) return `${name}: ${n} bonds`;
  }
  return void 0;
}

/** The CCD's standard InChI (`pdbx_chem_comp_descriptor`). */
export function ccdInchi(block: CifBlock): string | undefined {
  const descriptors = block.categories.pdbx_chem_comp_descriptor;
  if (!descriptors) return void 0;
  const type = descriptors.getField('type')!;
  const descriptor = descriptors.getField('descriptor')!;
  for (let i = 0; i < descriptors.rowCount; i++) {
    if (type.str(i) === 'InChI') return descriptor.str(i);
  }
  return void 0;
}
