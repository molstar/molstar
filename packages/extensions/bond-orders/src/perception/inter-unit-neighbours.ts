/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import type { Structure, Unit } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { InterUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';

/**
 * The covalent neighbours of each atom in other units of its model (a covalent ligand's bond to its protein), read
 * from the model as a structure (`IntraContext.modelStructure()`): a function of the model only, so perception stays
 * one of the model, the elements and the props (_spec §1.3). Keys and neighbours are model atoms (`ElementIndex`), not
 * unit indices: a neighbour is outside the unit being perceived, and the map is shared by every unit of the model.
 * Consumers look an atom up by `unit.elements[u]`. An atom can have several neighbours (one per altloc of its
 * partner, or a bridging atom). Each neighbour comes with its bond's order, as the model structure's built-in bonds
 * have it: perception never writes a bond between units, so the order is fixed (_spec §1.11).
 */
export type InterUnitNeighbours = ReadonlyMap<ElementIndex, readonly InterUnitNeighbour[]>;

export interface InterUnitNeighbour {
  /** the partner, a model atom */
  readonly element: ElementIndex;
  /** the bond's order */
  readonly order: number;
}

const NeighboursByModel = new WeakMap<Model, InterUnitNeighbours>();

export function getInterUnitNeighbours(model: Model, modelStructure: () => Structure): InterUnitNeighbours {
  const cached = NeighboursByModel.get(model);
  if (cached) return cached;

  const structure = modelStructure();
  const neighbours = new Map<ElementIndex, InterUnitNeighbour[]>();
  // a pair is listed under both units: each atom gets each partner once
  for (const e of sharedInterUnitBonds(structure, () => structure.interUnitBonds).edges) {
    if (!BondType.isCovalent(e.props.flag)) continue;
    const aI = structure.unitMap.get(e.unitA).elements[e.indexA];
    const neighbour = { element: structure.unitMap.get(e.unitB).elements[e.indexB], order: e.props.order };
    const list = neighbours.get(aI);
    if (list) list.push(neighbour);
    else neighbours.set(aI, [neighbour]);
  }
  NeighboursByModel.set(model, neighbours);
  return neighbours;
}

/** The built-in bonds between the units of a model's structure, with the units they were computed for. */
interface SharedBonds {
  readonly units: readonly Unit[];
  readonly bonds: InterUnitBonds;
}

const SharedByModel = new WeakMap<Model, SharedBonds>();

/**
 * The built-in bonds between the units of `structure`, computed once per model (_spec §2): the model structure that
 * perception reads (`getInterUnitNeighbours`) and the structure the viewer shows (the model, or an assembly of one
 * identity operator: same units, another operator object) share them. Reused only when `structure` has the stored
 * units (ids, elements, traits, operator matrices; the operator keys too when the model's bond table lists bonds per
 * operator) and the model structure's computation props (no `dynamicBonds`, no unit filters); otherwise `compute`.
 */
export function sharedInterUnitBonds(structure: Structure, compute: () => InterUnitBonds): InterUnitBonds {
  const model = singleModel(structure);
  if (!model || structure.dynamicBonds || structure.interBondsValidUnit || structure.interBondsValidUnitPair) {
    return compute();
  }
  const shared = SharedByModel.get(model);
  if (shared && haveSameUnits(model, shared.units, structure.units)) return shared.bonds;
  const bonds = compute();
  if (!shared) SharedByModel.set(model, { units: structure.units, bonds });
  return bonds;
}

function singleModel(structure: Structure): Model | undefined {
  const { models } = structure;
  return models.length === 1 ? models[0] : undefined;
}

function haveSameUnits(model: Model, a: readonly Unit[], b: readonly Unit[]) {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  // a bond table with bonds per operator matches them by operator key (inter-compute), not by matrix
  const byKey = !!IndexPairBonds.Provider.get(model)?.hasOperators;
  for (let i = 0, il = a.length; i < il; i++) {
    const u = a[i],
      v = b[i];
    if (u.id !== v.id || u.kind !== v.kind || u.traits !== v.traits) return false;
    const opU = u.conformation.operator,
      opV = v.conformation.operator;
    if (byKey && opU.key !== opV.key) return false;
    if (opU !== opV && !Mat4.areEqual(opU.matrix, opV.matrix, 0)) return false;
    if (!SortedArray.areEqual(u.elements, v.elements)) return false;
  }
  return true;
}
