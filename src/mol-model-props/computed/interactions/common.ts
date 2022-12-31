/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { InterUnitGraph } from '../../../mol-math/graph/inter-unit-graph';
import { Unit } from '../../../mol-model/structure';
import { AssignableArrayLike } from '../../../mol-util/type-helpers';
import { Features } from './features';
import { StructureElement } from '../../../mol-model/structure/structure';
import { IntMap } from '../../../mol-data/int';

export { InteractionsIntraContacts };
interface InteractionsIntraContacts extends IntAdjacencyGraph<Features.FeatureIndex, InteractionsIntraContacts.Props> {
    readonly elementsIndex: InteractionsIntraContacts.ElementsIndex
}
namespace InteractionsIntraContacts {
    export type Props = {
        readonly type: ArrayLike<InteractionType>
        readonly flag: AssignableArrayLike<InteractionFlag>
    }

    /** maps unit elements to contacts, range for unit element i is offsets[i] to offsets[i + 1] */
    export type ElementsIndex = {
        /** intra contact indices */
        readonly indices: ArrayLike<number>
        /** range for unit element i is offsets[i] to offsets[i + 1] */
        readonly offsets: ArrayLike<number>
    }

    /**
     * Note: assumes that feature members of a contact are non-overlapping
     */
    export function createElementsIndex(contacts: IntAdjacencyGraph<Features.FeatureIndex, Props>, features: Features, elementsCount: number) {
        const offsets = new Int32Array(elementsCount + 1);
        const bucketFill = new Int32Array(elementsCount);
        const bucketSizes = new Int32Array(elementsCount);
        const { members, offsets: featureOffsets } = features;

        for (let i = 0, il = contacts.edgeCount * 2; i < il; ++i) {
            const aI = contacts.a[i];
            const bI = contacts.b[i];
            if (aI > bI) continue;

            for (let j = featureOffsets[aI], jl = featureOffsets[aI + 1]; j < jl; ++j) {
                ++bucketSizes[members[j]];
            }
            for (let j = featureOffsets[bI], jl = featureOffsets[bI + 1]; j < jl; ++j) {
                ++bucketSizes[members[j]];
            }
        }

        let offset = 0;
        for (let i = 0; i < elementsCount; i++) {
            offsets[i] = offset;
            offset += bucketSizes[i];
        }
        offsets[elementsCount] = offset;

        const indices = new Int32Array(offset);
        for (let i = 0, il = contacts.edgeCount * 2; i < il; ++i) {
            const aI = contacts.a[i];
            const bI = contacts.b[i];
            if (aI > bI) continue;

            for (let j = featureOffsets[aI], jl = featureOffsets[aI + 1]; j < jl; ++j) {
                const m = members[j];
                const om = offsets[m] + bucketFill[m];
                indices[om] = i;
                ++bucketFill[m];
            }
            for (let j = featureOffsets[bI], jl = featureOffsets[bI + 1]; j < jl; ++j) {
                const m = members[j];
                const om = offsets[m] + bucketFill[m];
                indices[om] = i;
                ++bucketFill[m];
            }
        }

        return { indices, offsets };
    }
}

export { InteractionsInterContacts };
class InteractionsInterContacts extends InterUnitGraph<number, Features.FeatureIndex, InteractionsInterContacts.Props> {
    private readonly elementKeyIndex: Map<string, number[]>;

    getContactIndicesForElement(index: StructureElement.UnitIndex, unit: Unit): ReadonlyArray<number> {
        return this.elementKeyIndex.get(this.getElementKey(index, unit.id)) || [];
    }

    private getElementKey(index: StructureElement.UnitIndex, unitId: number): string {
        return `${index}|${unitId}`;
    }

    constructor(map: Map<number, InterUnitGraph.UnitPairEdges<number, Features.FeatureIndex, InteractionsInterContacts.Props>[]>, unitsFeatures: IntMap<Features>) {
        super(map);

        this.elementKeyIndex = new Map<string, number[]>();
        for (let i = 0, il = this.edges.length; i < il; ++i) {
            const { unitA, indexA } = this.edges[i];
            const { offsets, members } = unitsFeatures.get(unitA);
            for (let j = offsets[indexA], jl = offsets[indexA + 1]; j < jl; ++j) {
                const vertexKey = this.getElementKey(members[j], unitA);
                const e = this.elementKeyIndex.get(vertexKey);
                if (e === undefined) {
                    this.elementKeyIndex.set(vertexKey, [i]);
                } else {
                    e.push(i);
                }
            }
        }
    }
}
namespace InteractionsInterContacts {
    export type Props = { type: InteractionType, flag: InteractionFlag }
}

export enum InteractionFlag {
    None = 0,
    Filtered = 1,
}

export enum InteractionType {
    Unknown = 0,
    Ionic = 1,
    CationPi = 2,
    PiStacking = 3,
    HydrogenBond = 4,
    HalogenBond = 5,
    Hydrophobic = 6,
    MetalCoordination = 7,
    WeakHydrogenBond = 8,
}

export function interactionTypeLabel(type: InteractionType): string {
    switch (type) {
        case InteractionType.HydrogenBond:
            return 'Hydrogen Bond';
        case InteractionType.Hydrophobic:
            return 'Hydrophobic Contact';
        case InteractionType.HalogenBond:
            return 'Halogen Bond';
        case InteractionType.Ionic:
            return 'Ionic Interaction';
        case InteractionType.MetalCoordination:
            return 'Metal Coordination';
        case InteractionType.CationPi:
            return 'Cation-Pi Interaction';
        case InteractionType.PiStacking:
            return 'Pi Stacking';
        case InteractionType.WeakHydrogenBond:
            return 'Weak Hydrogen Bond';
        case InteractionType.Unknown:
            return 'Unknown Interaction';
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
    HydrophobicAtom = 8,
    WeakHydrogenDonor = 9,
    IonicTypePartner = 10,
    DativeBondPartner = 11,
    TransitionMetal = 12,
    IonicTypeMetal = 13
}

// to use with isolatedModules
export enum FeatureTypes {
    None = FeatureType.None,
    PositiveCharge = FeatureType.PositiveCharge,
    NegativeCharge = FeatureType.NegativeCharge,
    AromaticRing = FeatureType.AromaticRing,
    HydrogenDonor = FeatureType.HydrogenDonor,
    HydrogenAcceptor = FeatureType.HydrogenAcceptor,
    HalogenDonor = FeatureType.HalogenDonor,
    HalogenAcceptor = FeatureType.HalogenAcceptor,
    HydrophobicAtom = FeatureType.HydrophobicAtom,
    WeakHydrogenDonor = FeatureType.WeakHydrogenDonor,
    IonicTypePartner = FeatureType.IonicTypePartner,
    DativeBondPartner = FeatureType.DativeBondPartner,
    TransitionMetal = FeatureType.TransitionMetal,
    IonicTypeMetal = FeatureType.IonicTypeMetal
}

export function featureTypeLabel(type: FeatureType): string {
    switch (type) {
        case FeatureType.None:
            return 'None';
        case FeatureType.PositiveCharge:
            return 'Positive Charge';
        case FeatureType.NegativeCharge:
            return 'Negative Charge';
        case FeatureType.AromaticRing:
            return 'Aromatic Ring';
        case FeatureType.HydrogenDonor:
            return 'Hydrogen Donor';
        case FeatureType.HydrogenAcceptor:
            return 'Hydrogen Acceptor';
        case FeatureType.HalogenDonor:
            return 'Halogen Donor';
        case FeatureType.HalogenAcceptor:
            return 'Halogen Acceptor';
        case FeatureType.HydrophobicAtom:
            return 'HydrophobicAtom';
        case FeatureType.WeakHydrogenDonor:
            return 'Weak Hydrogen Donor';
        case FeatureType.IonicTypePartner:
            return 'Ionic Type Partner';
        case FeatureType.DativeBondPartner:
            return 'Dative Bond Partner';
        case FeatureType.TransitionMetal:
            return 'Transition Metal';
        case FeatureType.IonicTypeMetal:
            return 'Ionic Type Metal';
    }
}

export enum FeatureGroup {
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
            return 'None';
        case FeatureGroup.QuaternaryAmine:
            return 'Quaternary Amine';
        case FeatureGroup.TertiaryAmine:
            return 'Tertiary Amine';
        case FeatureGroup.Sulfonium:
            return 'Sulfonium';
        case FeatureGroup.SulfonicAcid:
            return 'Sulfonic Acid';
        case FeatureGroup.Sulfate:
            return 'Sulfate';
        case FeatureGroup.Phosphate:
            return 'Phosphate';
        case FeatureGroup.Halocarbon:
            return 'Halocarbon';
        case FeatureGroup.Guanidine:
            return 'Guanidine';
        case FeatureGroup.Acetamidine:
            return 'Acetamidine';
        case FeatureGroup.Carboxylate:
            return 'Carboxylate';
    }
}