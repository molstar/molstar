/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IntAdjacencyGraph } from '../../../mol-math/graph'
import { InterUnitGraph } from '../../../mol-math/graph/inter-unit-graph'
import { Unit } from '../../../mol-model/structure'

export type InteractionsIntraLinks = IntAdjacencyGraph<{ readonly type: ArrayLike<InteractionType> }>

export { InteractionsInterLinks }
type InteractionsInterLinks = InterUnitGraph<Unit, number, { type: InteractionType }>
namespace InteractionsInterLinks {
    export class Pair extends InterUnitGraph.UnitPairEdges<Unit, number, { type: InteractionType }> {}
    export type Info = InterUnitGraph.EdgeInfo<number, { type: InteractionType }>
}

export const enum InteractionType {
    Unknown = 0,
    IonicInteraction = 1,
    CationPi = 2,
    PiStacking = 3,
    HydrogenBond = 4,
    HalogenBond = 5,
    Hydrophobic = 6,
    MetalCoordination = 7,
    WeakHydrogenBond = 8,
    WaterHydrogenBond = 9,
    BackboneHydrogenBond = 10
}

export const enum FeatureType {
    None = 0,
    PositiveCharge = 1,
    NegativeCharge = 2,
    AromaticRing = 3,
    HydrogenDonor = 4,
    HydrogenAcceptor = 5,
    HalogenDonor = 6,
    HalogenAcceptor = 7,
    Hydrophobic = 8,
    WeakHydrogenDonor = 9,
    IonicTypePartner = 10,
    DativeBondPartner = 11,
    TransitionMetal = 12,
    IonicTypeMetal = 13
}

export const enum FeatureGroup {
    None = 0,
    QuaternaryAmine = 1,
    TertiaryAmine = 2,
    Sulfonium = 3,
    SulfonicAcid = 4,
    Sulfate = 5,
    Phosphate = 6,
    Halocarbon = 7,
    Guanidine = 8,
    Acetamidine = 9,
    Carboxylate = 10
}