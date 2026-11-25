/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import fs from 'node:fs';
import { CIF } from '@molstar/io/reader/cif';
import { parsePDB } from '@molstar/io/reader/pdb/parser';
import { trajectoryFromPDB } from '@molstar/model/formats/structure/pdb';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import { BondOrderPerceptionProvider } from '../provider.js';
import type { BondOrdersMode } from '../perception/perceive.js';

export function readFixture(name: string) {
  return fs.readFileSync(new URL(`./fixtures/${name}`, import.meta.url), 'utf8');
}

export async function modelFromPdb(text: string, id = 'test'): Promise<Model> {
  const parsed = await parsePDB(text, id).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromPDB(parsed.result).run()).representative;
}

/** A structure of the model with the perception provider in `mode`. */
export function perceived(model: Model, mode: BondOrdersMode = 'auto') {
  return Structure.ofModel(model, { bondProvider: BondOrderPerceptionProvider.create({ mode }) });
}

export function unitOf(s: Structure, compId: string) {
  const unit = s.units.find(
    (u) => Unit.isAtomic(u) && u.model.atomicHierarchy.atoms.label_comp_id.value(u.elements[0]) === compId,
  );
  if (!unit) throw new Error(`no unit of ${compId}`);
  return unit as Unit.Atomic;
}

export interface BondRow {
  a: string;
  b: string;
  order: number;
  aromatic: boolean;
}

/** The covalent bonds of a unit between named atoms, each once, `a < b`. */
export function bondRows(unit: Unit.Atomic, bonds: IntraUnitBonds = unit.bonds): BondRow[] {
  const { label_atom_id } = unit.model.atomicHierarchy.atoms;
  const { a, b, edgeProps } = bonds;
  const rows: BondRow[] = [];
  for (let i = 0; i < a.length; i++) {
    if (a[i] > b[i] || !BondType.isCovalent(edgeProps.flags[i])) continue;
    const na = label_atom_id.value(unit.elements[a[i]]);
    const nb = label_atom_id.value(unit.elements[b[i]]);
    rows.push({
      a: na < nb ? na : nb,
      b: na < nb ? nb : na,
      order: edgeProps.order[i],
      aromatic: BondType.is(edgeProps.flags[i], BondType.Flag.Aromatic),
    });
  }
  return rows;
}

/** The bonds of a CCD component definition (`chem_comp_bond`), as `BondRow`s. */
export async function ccdBondRows(text: string): Promise<BondRow[]> {
  const parsed = await CIF.parse(text).run();
  if (parsed.isError) throw new Error(parsed.message);
  const category = parsed.result.blocks[0].categories.chem_comp_bond;
  const id1 = category.getField('atom_id_1')!;
  const id2 = category.getField('atom_id_2')!;
  const valueOrder = category.getField('value_order')!;
  const aromaticFlag = category.getField('pdbx_aromatic_flag')!;
  const orders: Record<string, number> = { SING: 1, DOUB: 2, TRIP: 3, QUAD: 4 };
  const rows: BondRow[] = [];
  for (let i = 0; i < category.rowCount; i++) {
    const na = id1.str(i);
    const nb = id2.str(i);
    rows.push({
      a: na < nb ? na : nb,
      b: na < nb ? nb : na,
      order: orders[valueOrder.str(i).toUpperCase()] ?? 1,
      aromatic: aromaticFlag.str(i).toUpperCase() === 'Y',
    });
  }
  return rows;
}

/**
 * Perceived against expected, on the bonds the structure has: `missed` are expected multiple or aromatic bonds not
 * perceived as such (forward), `spurious` perceived multiple or aromatic bonds the expected ones do not have (reverse).
 * An aromatic bond matches on the flag alone (its Kekulé order may differ from the CCD's).
 */
export function compareBonds(perceived: BondRow[], expected: BondRow[]) {
  const key = (r: BondRow) => `${r.a}-${r.b}`;
  const expectedByKey = new Map(expected.map((r) => [key(r), r]));
  const perceivedByKey = new Map(perceived.map((r) => [key(r), r]));
  const matches = (p: BondRow, e: BondRow) => (e.aromatic ? p.aromatic : !p.aromatic && p.order === e.order);
  const missed: string[] = [];
  const spurious: string[] = [];
  for (const e of expected) {
    if (e.order === 1 && !e.aromatic) continue;
    const p = perceivedByKey.get(key(e));
    if (p && !matches(p, e)) missed.push(`${key(e)} expected ${e.aromatic ? 'aromatic' : e.order}, got ${p.order}`);
  }
  for (const p of perceived) {
    if (p.order === 1 && !p.aromatic) continue;
    const e = expectedByKey.get(key(p));
    if (!e || !matches(p, e)) spurious.push(`${key(p)} perceived ${p.aromatic ? 'aromatic' : p.order}`);
  }
  return { missed, spurious };
}
