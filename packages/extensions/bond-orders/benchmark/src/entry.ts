/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { BondOrderPerceptionProvider } from '@molstar/bond-orders-extension';
import { SortedArray } from '@molstar/core/data/int';
import { CIF } from '@molstar/io/reader/cif';
import type { CifBlock } from '@molstar/io/reader/cif/data-model';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { BondType } from '@molstar/model/model/structure/model/types';
import { bondRows, diagnostics, orderDiff } from './diagnostics.ts';
import { entryInput, isHydrogen } from './input.ts';
import type { Perceived } from './perceive.ts';
import { toMolBlock } from './sdf.ts';
import type { Tier } from './tiers.ts';

export interface EntryPerceived extends Perceived {
  /** the copy compared (plan §3.2) */
  copy: string;
  /** links to other residues, one cap each (plan §4.2) */
  caps: number;
}

/**
 * Tier C (plan §3.2): the component's first linked and complete copy in its entry, perceived in context
 * (`forceCompute`), beside the same copy with the entry's orders; both written with one cap per link (plan §4.2).
 */
export async function perceiveEntry(
  id: string,
  ccd: CifBlock,
  data: Uint8Array,
  tier: Tier,
  /** tier C: a linked copy only; a dataset's pairs (plan §6): the first complete copy, linked or not */
  linkedOnly = true,
): Promise<EntryPerceived | EntryExcluded> {
  const parsed = await CIF.parseBinary(data).run();
  if (parsed.isError) throw new Error(`${id}: ${parsed.message}`);
  const block = parsed.result.blocks[0];
  const trajectory = await trajectoryFromMmCIF(entryInput(block, { hydrogens: !tier.hydrogens, charges: true })).run();
  const model = (await trajectory.getFrameAtIndex(0)) as Model;
  const reference = Structure.ofModel(model);
  const copy = findCopy(id, ccd, reference, linkedOnly);
  if ('excluded' in copy) return copy;
  // after the copy: a remediated entry may no longer hold the component (and has no bonds for it)
  const dictionary = dictionaryDifference(id, ccd, block);
  if (dictionary) return { excluded: dictionary, linked: true };

  const perceived = Structure.ofModel(model, {
    bondProvider: BondOrderPerceptionProvider.create({ mode: 'forceCompute' }),
  });
  const isCap = (e: ElementIndex) => copy.caps.has(e);
  const subReference = substructure(reference, copy.elements);
  const subPerceived = substructure(perceived, copy.elements);
  const name = `${id} ${copy.label} tier ${tier.name}`;
  return {
    perceived: toMolBlock(subPerceived, name, isCap),
    reference: toMolBlock(subReference, `${id} ${copy.label} reference`, isCap),
    l0: orderDiff(bondRows(subReference, isCap), bondRows(subPerceived, isCap)),
    diagnostics: diagnostics(subPerceived, isCap),
    copy: copy.label,
    caps: copy.caps.size,
  };
}

/**
 * The file's bonds for the component against the CCD's, heavy atoms only, by atom pair and order (plan §3.2): a
 * component revised after the entry's release may differ.
 */
function dictionaryDifference(id: string, ccd: CifBlock, block: CifBlock): string | undefined {
  const heavy = new Set<string>();
  const atoms = ccd.categories.chem_comp_atom;
  for (let i = 0; i < atoms.rowCount; i++) {
    if (!isHydrogen(atoms.getField('type_symbol')!.str(i))) heavy.add(atoms.getField('atom_id')!.str(i));
  }
  const pairs = (category: CifBlock['categories'][string] | undefined, compId?: string) => {
    const set = new Set<string>();
    if (!category) return set;
    const comp = category.getField('comp_id');
    const a = category.getField('atom_id_1')!;
    const b = category.getField('atom_id_2')!;
    const order = category.getField('value_order')!;
    for (let i = 0; i < category.rowCount; i++) {
      if (compId && comp?.str(i) !== compId) continue;
      const [x, y] = [a.str(i), b.str(i)].sort();
      if (heavy.has(x) && heavy.has(y)) set.add(`${x}-${y} ${order.str(i).toUpperCase()}`);
    }
    return set;
  };
  const expected = pairs(ccd.categories.chem_comp_bond);
  const found = pairs(block.categories.chem_comp_bond, id);
  if (found.size === 0 && expected.size > 0) return 'entry has no chem_comp_bond for the component';
  const missing = [...expected].filter((p) => !found.has(p));
  const extra = [...found].filter((p) => !expected.has(p));
  if (missing.length || extra.length) {
    return `entry dictionary differs from the CCD (${[...missing.map((p) => `-${p}`), ...extra.map((p) => `+${p}`)].slice(0, 4).join(', ')})`;
  }
  return void 0;
}

/** `linked`: whether a copy covalently linked was found, its exclusion notwithstanding (the report's covalent ligands) */
export interface EntryExcluded {
  excluded: string;
  linked: boolean;
}

interface Copy {
  label: string;
  /** element indices per unit id: the copy's atoms and its caps */
  elements: Map<number, ElementIndex[]>;
  caps: Set<ElementIndex>;
}

/**
 * The first copy that is linked (a covalent bond to a heavy atom of another residue: metal coordination is not; any
 * copy when not `linkedOnly`) and complete (every heavy atom of the CCD but leaving atoms), with its link partners as
 * caps. A copy with a link that
 * replaces nothing of the definition is excluded (plan §3.2): the CCD's orders do not describe it.
 */
function findCopy(id: string, ccd: CifBlock, structure: Structure, linkedOnly: boolean): Copy | EntryExcluded {
  const required: string[] = [];
  const atoms = ccd.categories.chem_comp_atom;
  const leaving = atoms.getField('pdbx_leaving_atom_flag');
  for (let i = 0; i < atoms.rowCount; i++) {
    if (isHydrogen(atoms.getField('type_symbol')!.str(i))) continue;
    if (leaving?.str(i).toUpperCase() === 'Y') continue;
    required.push(atoms.getField('atom_id')!.str(i));
  }
  const substituents = substituentsByAtom(ccd);

  let found = 0;
  let linkedIncomplete = 0;
  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) continue;
    const { atomicHierarchy } = unit.model;
    const { label_seq_id, auth_seq_id } = atomicHierarchy.residues;
    const { label_atom_id, label_comp_id } = atomicHierarchy.atoms;
    const residueOf = atomicHierarchy.residueAtomSegments.index;
    const asymOf = atomicHierarchy.chains.label_asym_id;
    const chainOf = atomicHierarchy.chainAtomSegments.index;
    const { elements } = unit;
    for (let start = 0; start < elements.length; ) {
      const r = residueOf[elements[start]];
      let end = start;
      while (end < elements.length && residueOf[elements[end]] === r) end++;
      if (label_comp_id.value(elements[start]) === id) {
        found++;
        const names = new Set<string>();
        for (let i = start; i < end; i++) names.add(label_atom_id.value(elements[i]));
        const complete = required.every((n) => names.has(n));
        const caps = linkPartners(structure, unit, start, end);
        if (caps.length && !complete) linkedIncomplete++;
        const unreplaced = caps.length && complete ? unreplacedLink(caps, names, substituents, label_atom_id) : void 0;
        if (unreplaced) {
          return { excluded: `link replaces no H or leaving atom of the CCD (${unreplaced})`, linked: true };
        }
        if (complete && (caps.length || !linkedOnly)) {
          const copyElements = new Map<number, ElementIndex[]>();
          copyElements.set(
            unit.id,
            Array.from({ length: end - start }, (_, k) => elements[start + k]),
          );
          const capSet = new Set<ElementIndex>();
          for (const { unitId, e } of caps) {
            // one cap per partner, even bonded to two atoms of the copy
            if (capSet.has(e)) continue;
            capSet.add(e);
            const list = copyElements.get(unitId) ?? [];
            list.push(e);
            copyElements.set(unitId, list);
          }
          const label = `${asymOf.value(chainOf[elements[start]])} ${label_seq_id.value(r)} (auth ${auth_seq_id.value(r)})`;
          return { label, elements: copyElements, caps: capSet };
        }
      }
      start = end;
    }
  }
  if (!found) return { excluded: 'component not found in its entry', linked: false };
  if (!linkedOnly) return { excluded: 'entry incomplete', linked: linkedIncomplete > 0 };
  if (linkedIncomplete) return { excluded: 'representative entry incomplete', linked: true };
  return { excluded: 'not linked in its entry', linked: false };
}

/**
 * Per heavy atom of the CCD definition, how many links it can take without changing the definition's orders: its H
 * (`substitutable`), and its leaving heavy neighbours (`leaving`), which count only when absent from the copy.
 */
function substituentsByAtom(ccd: CifBlock) {
  const atoms = ccd.categories.chem_comp_atom;
  const hydrogen = new Set<string>();
  const leaving = new Set<string>();
  const leavingFlag = atoms.getField('pdbx_leaving_atom_flag');
  for (let i = 0; i < atoms.rowCount; i++) {
    const name = atoms.getField('atom_id')!.str(i);
    if (isHydrogen(atoms.getField('type_symbol')!.str(i))) hydrogen.add(name);
    else if (leavingFlag?.str(i).toUpperCase() === 'Y') leaving.add(name);
  }
  const byAtom = new Map<string, { substitutable: number; leaving: string[] }>();
  const get = (name: string) => {
    let s = byAtom.get(name);
    if (!s) byAtom.set(name, (s = { substitutable: 0, leaving: [] }));
    return s;
  };
  const bonds = ccd.categories.chem_comp_bond;
  if (!bonds) return byAtom;
  const a = bonds.getField('atom_id_1')!;
  const b = bonds.getField('atom_id_2')!;
  for (let i = 0; i < bonds.rowCount; i++) {
    for (const [x, y] of [
      [a.str(i), b.str(i)],
      [b.str(i), a.str(i)],
    ]) {
      if (hydrogen.has(x)) continue;
      if (hydrogen.has(y)) get(x).substitutable++;
      else if (leaving.has(y)) get(x).leaving.push(y);
    }
  }
  return byAtom;
}

/**
 * An atom of the copy whose links outnumber what they can replace, its H in the CCD and its leaving neighbours missing
 * from the copy, if any. Condensations (peptide, glycosidic and phosphodiester bonds) replace one or the other. Mostly
 * an addition at an atom without H (a boronate, a carbonyl C): the CCD has the unbound form; otherwise a spurious link
 * or a leaving atom the CCD does not flag. An addition at an atom with H (an acrylamide's CH2) passes.
 */
function unreplacedLink(
  links: LinkPartner[],
  names: Set<string>,
  substituents: ReturnType<typeof substituentsByAtom>,
  label_atom_id: { value(e: ElementIndex): string },
) {
  const count = new Map<string, number>();
  for (const { inner } of links) {
    const name = label_atom_id.value(inner);
    count.set(name, (count.get(name) ?? 0) + 1);
  }
  for (const [name, n] of count) {
    const s = substituents.get(name);
    const replaceable = (s?.substitutable ?? 0) + (s?.leaving.filter((l) => !names.has(l)).length ?? 0);
    if (n > replaceable) return name;
  }
  return void 0;
}

interface LinkPartner {
  unitId: number;
  /** the partner, a heavy atom of another residue */
  e: ElementIndex;
  /** the copy's atom it is bonded to */
  inner: ElementIndex;
}

/** The covalent links of the residue at unit indices `start..end` to heavy atoms of other residues. */
function linkPartners(structure: Structure, unit: Unit.Atomic, start: number, end: number): LinkPartner[] {
  const partners = new Map<string, LinkPartner>();
  const { type_symbol } = unit.model.atomicHierarchy.atoms;
  const add = (inner: ElementIndex, u: Unit, e: ElementIndex, flags: number) => {
    if (!BondType.isCovalent(flags) || isHydrogen(type_symbol.value(e))) return;
    partners.set(`${inner}-${u.id}:${e}`, { unitId: u.id, e, inner });
  };
  const { offset, b, edgeProps } = unit.bonds;
  for (let i = start; i < end; i++) {
    for (let t = offset[i]; t < offset[i + 1]; t++) {
      if (b[t] < start || b[t] >= end) add(unit.elements[i], unit, unit.elements[b[t]], edgeProps.flags[t]);
    }
  }
  for (const edge of structure.interUnitBonds.edges) {
    if (edge.unitA !== unit.id || edge.indexA < start || edge.indexA >= end) continue;
    const other = structure.unitMap.get(edge.unitB);
    add(unit.elements[edge.indexA], other, other.elements[edge.indexB], edge.props.flag);
  }
  return [...partners.values()];
}

/** The structure restricted to `elements` (per unit id), its bonds the parent's (`childBonds: 'root'`). */
function substructure(structure: Structure, elements: Map<number, ElementIndex[]>): Structure {
  const units: Unit[] = [];
  for (const [unitId, list] of elements) {
    const unit = structure.unitMap.get(unitId);
    units.push(unit.getChild(SortedArray.ofUnsortedArray(list)));
  }
  return Structure.create(units, { parent: structure });
}
