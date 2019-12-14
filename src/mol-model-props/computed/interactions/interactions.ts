/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { addUnitHydrogenDonors, addUnitWeakHydrogenDonors, addUnitHydrogenAcceptors, addHydrogenBonds, HydrogenBondsParams } from './hydrogen-bonds';
import { Features, FeaturesBuilder } from './features';
import { ValenceModelProvider } from '../valence-model';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { RuntimeContext } from '../../../mol-task';

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

export { InteractionsBuilder }

interface InteractionsBuilder {
    add: (indexA: number, indexB: number, type: InteractionType) => void
    getInteractions: () => Interactions
}

namespace InteractionsBuilder {
    export function create(features: Features, elementsCount: number): InteractionsBuilder {
        const aIndices: number[] = []
        const bIndices: number[] = []
        const types: number[] = []

        return {
            add(indexA: number, indexB: number, type: InteractionType) {
                aIndices[aIndices.length] = indexA
                bIndices[bIndices.length] = indexB
                types[types.length] = type
            },
            getInteractions() {
                const builder = new IntAdjacencyGraph.EdgeBuilder(features.count, aIndices, bIndices)
                const _types = new Int8Array(builder.slotCount) as ArrayLike<InteractionType>
                for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
                    builder.addNextEdge()
                    builder.assignProperty(_types, types[i])
                }
                const links = builder.createGraph({ types: _types })

                const elementsIndex = Features.createElementsIndex(features, elementsCount)

                return {
                    links,
                    features,
                    elementsIndex,
                    getLinkIndex: (indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex) => {
                        // TODO quadratic runtime... when both indices are part of more than one feature
                        const { indices, offsets } = elementsIndex
                        for (let i = offsets[indexA], il = offsets[indexA + 1]; i < il; ++i) {
                            const fA = indices[i]
                            for (let j = offsets[indexB], jl = offsets[indexB + 1]; j < jl; ++j) {
                                const fB = indices[j]
                                let l = links.getDirectedEdgeIndex(fA, fB)
                                if (l !== -1) return l
                            }
                        }
                        return -1
                    }
                }
            }
        }
    }
}

export type InteractionsLinks = IntAdjacencyGraph<{ readonly types: ArrayLike<InteractionType> }>

export interface Interactions {
    links: InteractionsLinks
    features: Features
    elementsIndex: Features.ElementsIndex
    getLinkIndex: (indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex) => number
}

export const InteractionsParams = {
    ...HydrogenBondsParams
}
export type InteractionsParams = typeof InteractionsParams
export type InteractionsProps = PD.Values<InteractionsParams>

export async function calcInteractions(runtime: RuntimeContext, structure: Structure, props: Partial<InteractionsProps>) {
    const p = { ...PD.getDefaultValues(InteractionsParams), ...props }
    await ValenceModelProvider.attach(structure).runInContext(runtime)
    const map = new Map<number, Interactions>()
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const u = structure.units[i]
        if (Unit.isAtomic(u)) {
            const interactions = calcIntraUnitInteractions(structure, u, p)
            map.set(u.id, interactions)
        }
    }
    return map
}

function calcIntraUnitInteractions(structure: Structure, unit: Unit.Atomic, props: InteractionsProps) {

    const featuresBuilder = FeaturesBuilder.create()
    addUnitHydrogenDonors(structure, unit, featuresBuilder)
    addUnitWeakHydrogenDonors(structure, unit, featuresBuilder)
    addUnitHydrogenAcceptors(structure, unit, featuresBuilder)

    const features = featuresBuilder.getFeatures()

    const interactionsBuilder = InteractionsBuilder.create(features, unit.elements.length)
    addHydrogenBonds(structure, unit, features, interactionsBuilder, props)

    return interactionsBuilder.getInteractions()
}