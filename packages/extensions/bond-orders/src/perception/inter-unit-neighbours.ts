/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import type { ElementIndex } from '@molstar/model/model/structure/model/indexing';
import { BondType } from '@molstar/model/model/structure/model/types';

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
  for (const e of structure.interUnitBonds.edges) {
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
