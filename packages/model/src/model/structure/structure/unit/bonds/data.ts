/**
 * Copyright (c) 2017-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { BondType } from '../../../model/types.js';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { Unit } from '../../unit.js';
import type { StructureElement } from '../../element.js';
import type { Bond } from '../bonds.js';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';

export type IntraUnitBondProps = {
  /** Can remap even with `dynamicBonds` on, e.g., for water molecules */
  readonly canRemap?: boolean;
  /** Can be cached in `ElementSetIntraBondCache` */
  readonly cacheable?: boolean;
};

export type IntraUnitBonds = IntAdjacencyGraph<
  StructureElement.UnitIndex,
  {
    readonly order: ArrayLike<number>;
    readonly flags: ArrayLike<BondType.Flag>;
    readonly key: ArrayLike<number>;
  },
  IntraUnitBondProps
>;

export namespace IntraUnitBonds {
  export const Empty: IntraUnitBonds = IntAdjacencyGraph.create([], [], [], 0, { flags: [], order: [], key: [] });

  /**
   * The bonds of `atomCount` atoms from a list of pairs, each once (`atomA[i]`–`atomB[i]`, unit indices), with their
   * order, flags and key (-1: none): what the built-in computation builds, and what a bond provider can use.
   */
  export function fromPairs(
    atomA: ArrayLike<StructureElement.UnitIndex>,
    atomB: ArrayLike<StructureElement.UnitIndex>,
    order: ArrayLike<number>,
    flags: ArrayLike<number>,
    key: ArrayLike<number>,
    atomCount: number,
    props?: IntraUnitBondProps,
  ): IntraUnitBonds {
    const builder = new IntAdjacencyGraph.EdgeBuilder(atomCount, atomA, atomB);
    const flagsOut = new Uint16Array(builder.slotCount);
    const orderOut = new Int8Array(builder.slotCount);
    const keyOut = new Uint32Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
      builder.addNextEdge();
      builder.assignProperty(flagsOut, flags[i]);
      builder.assignProperty(orderOut, order[i]);
      builder.assignProperty(keyOut, key[i]);
    }
    return builder.createGraph({ flags: flagsOut, order: orderOut, key: keyOut }, props);
  }
}

export type InterUnitEdgeProps = { readonly order: number; readonly flag: BondType.Flag; readonly key: number };

export class InterUnitBonds extends InterUnitGraph<number, StructureElement.UnitIndex, InterUnitEdgeProps> {
  /** Get inter-unit bond given a bond-location */
  getBondFromLocation(l: Bond.Location) {
    return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit)
      ? this.getEdge(l.aIndex, l.aUnit.id, l.bIndex, l.bUnit.id)
      : undefined;
  }

  /** Get inter-unit bond index given a bond-location */
  getBondIndexFromLocation(l: Bond.Location) {
    return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit)
      ? this.getEdgeIndex(l.aIndex, l.aUnit.id, l.bIndex, l.bUnit.id)
      : -1;
  }
}

export namespace InterUnitBonds {
  export class UnitPairBonds extends InterUnitGraph.UnitPairEdges<
    number,
    StructureElement.UnitIndex,
    InterUnitEdgeProps
  > {}
  export type BondInfo = InterUnitGraph.EdgeInfo<StructureElement.UnitIndex, InterUnitEdgeProps>;
}
