/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit } from '../../../mol-model/structure';
import { RuntimeContext } from '../../../mol-task';
import { Features, FeaturesBuilder } from './features';
import { ValenceModelProvider } from '../valence-model';
import { InteractionsIntraContacts, InteractionsInterContacts, FeatureType } from './common';
import { IntraContactsBuilder, InterContactsBuilder } from './contacts-builder';
import { IntMap } from '../../../mol-data/int';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { addUnitContacts, ContactTester, addStructureContacts, ContactProvider, ContactsParams, ContactsProps } from './contacts';
import { HalogenDonorProvider, HalogenAcceptorProvider, HalogenBondsProvider } from './halogen-bonds';
import { HydrogenDonorProvider, WeakHydrogenDonorProvider, HydrogenAcceptorProvider, HydrogenBondsProvider, WeakHydrogenBondsProvider } from './hydrogen-bonds';
import { NegativChargeProvider, PositiveChargeProvider, AromaticRingProvider, IonicProvider, PiStackingProvider, CationPiProvider } from './charged';
import { HydrophobicAtomProvider, HydrophobicProvider } from './hydrophobic';
import { SetUtils } from '../../../mol-util/set';
import { MetalCoordinationProvider, MetalProvider, MetalBindingProvider } from './metal';
import { refineInteractions } from './refine';

export { Interactions }

interface Interactions {
    /** Features of each unit */
    unitsFeatures: IntMap<Features>
    /** Interactions of each unit */
    unitsContacts: IntMap<InteractionsIntraContacts>
    /** Interactions between units */
    contacts: InteractionsInterContacts
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
}

const FeatureProviders = [
    HydrogenDonorProvider, WeakHydrogenDonorProvider, HydrogenAcceptorProvider,
    NegativChargeProvider, PositiveChargeProvider, AromaticRingProvider,
    HalogenDonorProvider, HalogenAcceptorProvider,
    HydrophobicAtomProvider,
    MetalProvider, MetalBindingProvider,
]

const LinkProviders = {
    'ionic': IonicProvider,
    'pi-stacking': PiStackingProvider,
    'cation-pi': CationPiProvider,
    'halogen-bonds': HalogenBondsProvider,
    'hydrogen-bonds': HydrogenBondsProvider,
    'weak-hydrogen-bonds': WeakHydrogenBondsProvider,
    'hydrophobic': HydrophobicProvider,
    'metal-coordination': MetalCoordinationProvider,
}
type LinkProviders = typeof LinkProviders

function getProvidersParams() {
    const params: { [k in keyof LinkProviders]: PD.Group<LinkProviders[k]['params']> } = Object.create(null)
    Object.keys(LinkProviders).forEach(k => {
        (params as any)[k] = PD.Group(LinkProviders[k as keyof LinkProviders].params)
    })
    return params
}
export const InteractionsParams = {
    types: PD.MultiSelect([
        'ionic',
        'cation-pi',
        'pi-stacking',
        'hydrogen-bonds',
        'halogen-bonds',
        // 'hydrophobic',
        'metal-coordination',
        // 'weak-hydrogen-bonds',
    ], PD.objectToOptions(LinkProviders)),
    contacts: PD.Group(ContactsParams, { isFlat: true }),
    ...getProvidersParams()
}
export type InteractionsParams = typeof InteractionsParams
export type InteractionsProps = PD.Values<InteractionsParams>

export async function computeInteractions(runtime: RuntimeContext, structure: Structure, props: Partial<InteractionsProps>): Promise<Interactions> {
    const p = { ...PD.getDefaultValues(InteractionsParams), ...props }
    await ValenceModelProvider.attach(structure).runInContext(runtime)

    const linkProviders: ContactProvider<any>[] = []
    Object.keys(LinkProviders).forEach(k => {
        if (p.types.includes(k)) linkProviders.push(LinkProviders[k as keyof typeof LinkProviders])
    })
    const linkTesters = linkProviders.map(l => l.createTester(p[l.name as keyof InteractionsProps]))

    const requiredFeatures = new Set<FeatureType>()
    linkTesters.forEach(l => SetUtils.add(requiredFeatures, l.requiredFeatures))
    const featureProviders = FeatureProviders.filter(f => SetUtils.areIntersecting(requiredFeatures, f.types))

    const unitsFeatures = IntMap.Mutable<Features>()
    const unitsContacts = IntMap.Mutable<InteractionsIntraContacts>()

    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const group = structure.unitSymmetryGroups[i]
        if (runtime.shouldUpdate) {
            await runtime.update({ message: 'computing interactions', current: i, max: il })
        }
        const features = findUnitFeatures(structure, group.units[0], featureProviders)
        const intraUnitLinks = findIntraUnitContacts(structure, group.units[0], features, linkTesters, p.contacts)
        for (let j = 0, jl = group.units.length; j < jl; ++j) {
            const u = group.units[j]
            unitsFeatures.set(u.id, features)
            unitsContacts.set(u.id, intraUnitLinks)
        }
    }

    const contacts = findInterUnitContacts(structure, unitsFeatures, linkTesters, p.contacts)

    const interactions = { unitsFeatures, unitsContacts, contacts }
    refineInteractions(structure, interactions)
    return interactions
}

function findUnitFeatures(structure: Structure, unit: Unit, featureProviders: Features.Provider[]) {
    const count = unit.elements.length
    const featuresBuilder = FeaturesBuilder.create(count, count / 2)
    if (Unit.isAtomic(unit)) {
        for (const fp of featureProviders) {
            fp.add(structure, unit, featuresBuilder)
        }
    }
    return featuresBuilder.getFeatures(count)
}

function findIntraUnitContacts(structure: Structure, unit: Unit, features: Features, contactTesters: ReadonlyArray<ContactTester>, props: ContactsProps) {
    const linksBuilder = IntraContactsBuilder.create(features, unit.elements.length)
    if (Unit.isAtomic(unit)) {
        addUnitContacts(structure, unit, features, linksBuilder, contactTesters, props)
    }
    return linksBuilder.getContacts()
}

function findInterUnitContacts(structure: Structure, unitsFeatures: IntMap<Features>, contactTesters: ReadonlyArray<ContactTester>, props: ContactsProps) {
    const builder = InterContactsBuilder.create()

    const maxDistance = Math.max(...contactTesters.map(t => t.maxDistance))

    const lookup = structure.lookup3d;
    const imageCenter = Vec3.zero();

    for (const unitA of structure.units) {
        if (!Unit.isAtomic(unitA)) continue;

        const featuresA = unitsFeatures.get(unitA.id)

        const bs = unitA.lookup3d.boundary.sphere;
        Vec3.transformMat4(imageCenter, bs.center, unitA.conformation.operator.matrix);
        const closeUnits = lookup.findUnitIndices(imageCenter[0], imageCenter[1], imageCenter[2], bs.radius + maxDistance);

        for (let i = 0; i < closeUnits.count; i++) {
            const unitB = structure.units[closeUnits.indices[i]];
            if (!Unit.isAtomic(unitB) || unitA.id >= unitB.id || !Structure.validUnitPair(structure, unitA, unitB)) continue;

            const featuresB = unitsFeatures.get(unitB.id)

            if (unitB.elements.length >= unitA.elements.length) {
                addStructureContacts(structure, unitA, featuresA, unitB, featuresB, builder, contactTesters, props)
            } else {
                addStructureContacts(structure, unitB, featuresB, unitA, featuresA, builder, contactTesters, props)
            }
        }
    }

    return builder.getContacts()
}