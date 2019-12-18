/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit } from '../../../mol-model/structure';
import { RuntimeContext } from '../../../mol-task';
import { addUnitHydrogenDonors, addUnitWeakHydrogenDonors, addUnitHydrogenAcceptors, addUnitHydrogenBonds, HydrogenBondsParams, addStructureHydrogenBonds } from './hydrogen-bonds';
import { Features, FeaturesBuilder } from './features';
import { ValenceModelProvider } from '../valence-model';
import { InteractionsIntraLinks, InteractionsInterLinks, InteractionType } from './common';
import { IntraLinksBuilder, InterLinksBuilder } from './builder';
import { IntMap } from '../../../mol-data/int';
import { Vec3 } from '../../../mol-math/linear-algebra';

export { Interactions }

interface Interactions {
    /** Features of each unit */
    unitsFeatures: IntMap<Features>
    /** Interactions of each unit */
    unitsLinks: IntMap<InteractionsIntraLinks>
    /** Interactions between units */
    links: InteractionsInterLinks
}

namespace Interactions {
    export interface Location {
        readonly kind: 'interaction-location'
        interactions: Interactions
        unitA: Unit
        /** Index into features of unitA */
        indexA: number
        unitB: Unit
        /** Index into features of unitB */
        indexB: number
    }

    export function Location(interactions?: Interactions, unitA?: Unit, indexA?: number, unitB?: Unit, indexB?: number): Location {
        return { kind: 'interaction-location', interactions: interactions as any, unitA: unitA as any, indexA: indexA as any, unitB: unitB as any, indexB: indexB as any };
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'interaction-location';
    }

    export function areLocationsEqual(locA: Location, locB: Location) {
        return (
            locA.interactions === locB.interactions &&
            locA.indexA === locB.indexA && locA.indexB === locB.indexB &&
            locA.unitA === locB.unitA && locA.unitB === locB.unitB
        )
    }

    export interface Loci {
        readonly kind: 'interaction-loci'
        readonly structure: Structure
        readonly interactions: Interactions
        readonly links: ReadonlyArray<{
            unitA: Unit
            /** Index into features of unitA */
            indexA: number
            unitB: Unit
            /** Index into features of unitB */
            indexB: number
        }>
    }

    export function Loci(structure: Structure, interactions: Interactions, links: Loci['links']): Loci {
        return { kind: 'interaction-loci', structure, interactions, links };
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'interaction-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        if (a.structure !== b.structure) return false
        if (a.interactions !== b.interactions) return false
        if (a.links.length !== b.links.length) return false
        for (let i = 0, il = a.links.length; i < il; ++i) {
            const linkA = a.links[i]
            const linkB = b.links[i]
            if (linkA.unitA !== linkB.unitA) return false
            if (linkA.unitB !== linkB.unitB) return false
            if (linkA.indexA !== linkB.indexA) return false
            if (linkA.indexB !== linkB.indexB) return false
        }
        return true
    }

    export function isLociEmpty(loci: Loci) {
        return loci.links.length === 0 ? true : false
    }

    export function typeLabel(type: InteractionType): string {
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
}

export const InteractionsParams = {
    hydrogenBonds: PD.Group(HydrogenBondsParams),
}
export type InteractionsParams = typeof InteractionsParams
export type InteractionsProps = PD.Values<InteractionsParams>

export async function computeInteractions(runtime: RuntimeContext, structure: Structure, props: Partial<InteractionsProps>) {
    const p = { ...PD.getDefaultValues(InteractionsParams), ...props }
    await ValenceModelProvider.attach(structure).runInContext(runtime)

    const unitsFeatures = IntMap.Mutable<Features>()
    const unitsLinks = IntMap.Mutable<InteractionsIntraLinks>()

    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const group = structure.unitSymmetryGroups[i]
        const d = findIntraUnitLinksAndFeatures(structure, group.units[0], p)
        for (let j = 0, jl = group.units.length; j < jl; ++j) {
            const u = group.units[j]
            unitsFeatures.set(u.id, d.features)
            unitsLinks.set(u.id, d.links)
        }
    }

    const links = findInterUnitLinks(structure, unitsFeatures, p)

    return { unitsFeatures, unitsLinks, links }
}

function findIntraUnitLinksAndFeatures(structure: Structure, unit: Unit, props: InteractionsProps) {

    const featuresBuilder = FeaturesBuilder.create()
    if (Unit.isAtomic(unit)) {
        addUnitHydrogenDonors(structure, unit, featuresBuilder)
        addUnitWeakHydrogenDonors(structure, unit, featuresBuilder)
        addUnitHydrogenAcceptors(structure, unit, featuresBuilder)
    }
    const features = featuresBuilder.getFeatures()

    const linksBuilder = IntraLinksBuilder.create(features, unit.elements.length)
    if (Unit.isAtomic(unit)) {
        addUnitHydrogenBonds(structure, unit, features, linksBuilder, props.hydrogenBonds)
    }

    return { features, links: linksBuilder.getLinks() }
}

const MAX_RADIUS = 5

function findInterUnitLinks(structure: Structure, unitsFeatures: IntMap<Features>, props: InteractionsProps) {
    const builder = InterLinksBuilder.create()

    const lookup = structure.lookup3d;
    const imageCenter = Vec3.zero();

    for (const unitA of structure.units) {
        if (!Unit.isAtomic(unitA)) continue;

        const featuresA = unitsFeatures.get(unitA.id)

        const bs = unitA.lookup3d.boundary.sphere;
        Vec3.transformMat4(imageCenter, bs.center, unitA.conformation.operator.matrix);
        const closeUnits = lookup.findUnitIndices(imageCenter[0], imageCenter[1], imageCenter[2], bs.radius + MAX_RADIUS);

        for (let i = 0; i < closeUnits.count; i++) {
            const unitB = structure.units[closeUnits.indices[i]];
            if (!Unit.isAtomic(unitB) || unitA.id >= unitB.id || !Structure.validUnitPair(structure, unitA, unitB)) continue;

            const featuresB = unitsFeatures.get(unitB.id)

            if (unitB.elements.length >= unitA.elements.length) {
                addStructureHydrogenBonds(structure, unitA, featuresA, unitB, featuresB, builder, props.hydrogenBonds)
            } else {
                addStructureHydrogenBonds(structure, unitB, featuresB, unitA, featuresA, builder, props.hydrogenBonds)
            }
        }
    }

    return builder.getLinks()
}