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
 * partner, or a bridging atom).
 */
export type InterUnitNeighbours = ReadonlyMap<ElementIndex, readonly ElementIndex[]>;

const NeighboursByModel = new WeakMap<Model, InterUnitNeighbours>();

export function getInterUnitNeighbours(model: Model, modelStructure: () => Structure): InterUnitNeighbours {
  const cached = NeighboursByModel.get(model);
  if (cached) return cached;

  const structure = modelStructure();
  const neighbours = new Map<ElementIndex, ElementIndex[]>();
  // a pair is listed under both units: each atom gets each partner once
  for (const e of structure.interUnitBonds.edges) {
    if (!BondType.isCovalent(e.props.flag)) continue;
    const aI = structure.unitMap.get(e.unitA).elements[e.indexA];
    const bI = structure.unitMap.get(e.unitB).elements[e.indexB];
    const list = neighbours.get(aI);
    if (list) list.push(bI);
    else neighbours.set(aI, [bI]);
  }
  NeighboursByModel.set(model, neighbours);
  return neighbours;
}
