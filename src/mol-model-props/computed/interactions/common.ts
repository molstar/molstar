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

export function interactionTypeLabel(type: InteractionType): string {
    switch (type) {
        case InteractionType.HydrogenBond:
        case InteractionType.WaterHydrogenBond:
        case InteractionType.BackboneHydrogenBond:
            return 'Hydrogen Bond'
        case InteractionType.Hydrophobic:
            return 'Hydrophobic Contact'
        case InteractionType.HalogenBond:
            return 'Halogen Bond'
        case InteractionType.IonicInteraction:
            return 'Ionic Interaction'
        case InteractionType.MetalCoordination:
            return 'Metal Coordination'
        case InteractionType.CationPi:
            return 'Cation-Pi Interaction'
        case InteractionType.PiStacking:
            return 'Pi Stacking'
        case InteractionType.WeakHydrogenBond:
            return 'Weak Hydrogen Bond'
        case InteractionType.Unknown:
            return 'Unknown Interaction'
    }
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

export function featureTypeLabel(type: FeatureType): string {
    switch (type) {
        case FeatureType.None:
            return 'None'
        case FeatureType.PositiveCharge:
            return 'Positive Charge'
        case FeatureType.NegativeCharge:
            return 'Negative Charge'
        case FeatureType.AromaticRing:
            return 'Aromatic Ring'
        case FeatureType.HydrogenDonor:
            return 'Hydrogen Donor'
        case FeatureType.HydrogenAcceptor:
            return 'Hydrogen Acceptor'
        case FeatureType.HalogenDonor:
            return 'Halogen Donor'
        case FeatureType.HalogenAcceptor:
            return 'Halogen Acceptor'
        case FeatureType.Hydrophobic:
            return 'Hydrophobic'
        case FeatureType.WeakHydrogenDonor:
            return 'Weak Hydrogen Donor'
        case FeatureType.IonicTypePartner:
            return 'Ionic Type Partner'
        case FeatureType.DativeBondPartner:
            return 'Dative Bond Partner'
        case FeatureType.TransitionMetal:
            return 'Transition Metal'
        case FeatureType.IonicTypeMetal:
            return 'Ionic Type Metal'
    }
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

export function featureGroupLabel(group: FeatureGroup): string {
    switch (group) {
        case FeatureGroup.None:
            return 'None'
        case FeatureGroup.QuaternaryAmine:
            return 'Quaternary Amine'
        case FeatureGroup.TertiaryAmine:
            return 'Tertiary Amine'
        case FeatureGroup.Sulfonium:
            return 'Sulfonium'
        case FeatureGroup.SulfonicAcid:
            return 'Sulfonic Acid'
        case FeatureGroup.Sulfate:
            return 'Sulfate'
        case FeatureGroup.Phosphate:
            return 'Phosphate'
        case FeatureGroup.Halocarbon:
            return 'Halocarbon'
        case FeatureGroup.Guanidine:
            return 'Guanidine'
        case FeatureGroup.Acetamidine:
            return 'Acetamidine'
        case FeatureGroup.Carboxylate:
            return 'Carboxylate'
    }
}