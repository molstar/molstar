/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, Bond } from '../../../mol-model/structure';
import { Features, FeaturesBuilder } from './features';
import { ValenceModelProvider } from '../valence-model';
import { InteractionsIntraContacts, InteractionsInterContacts, FeatureType, interactionTypeLabel } from './common';
import { IntraContactsBuilder, InterContactsBuilder } from './contacts-builder';
import { IntMap } from '../../../mol-data/int';
import { addUnitContacts, ContactTester, addStructureContacts, ContactsParams, ContactsProps } from './contacts';
import { HalogenDonorProvider, HalogenAcceptorProvider, HalogenBondsProvider } from './halogen-bonds';
import { HydrogenDonorProvider, WeakHydrogenDonorProvider, HydrogenAcceptorProvider, HydrogenBondsProvider, WeakHydrogenBondsProvider } from './hydrogen-bonds';
import { NegativChargeProvider, PositiveChargeProvider, AromaticRingProvider, IonicProvider, PiStackingProvider, CationPiProvider } from './charged';
import { HydrophobicAtomProvider, HydrophobicProvider } from './hydrophobic';
import { SetUtils } from '../../../mol-util/set';
import { MetalCoordinationProvider, MetalProvider, MetalBindingProvider } from './metal';
import { refineInteractions } from './refine';
import { CustomProperty } from '../../common/custom-property';
import { DataLocation } from '../../../mol-model/location';
import { CentroidHelper } from '../../../mol-math/geometry/centroid-helper';
import { Sphere3D } from '../../../mol-math/geometry';
import { DataLoci } from '../../../mol-model/loci';
import { bondLabel, LabelGranularity } from '../../../mol-theme/label';
import { ObjectKeys } from '../../../mol-util/type-helpers';

export { Interactions };

interface Interactions {
    /** Features of each unit */
    unitsFeatures: IntMap<Features>
    /** Interactions of each unit */
    unitsContacts: IntMap<InteractionsIntraContacts>
    /** Interactions between units */
    contacts: InteractionsInterContacts
}

namespace Interactions {
    type StructureInteractions = { readonly structure: Structure, readonly interactions: Interactions }

    export interface Element {
        unitA: Unit
        /** Index into features of unitA */
        indexA: Features.FeatureIndex
        unitB: Unit
        /** Index into features of unitB */
        indexB: Features.FeatureIndex
    }

    export interface Location extends DataLocation<StructureInteractions, Element> {}

    export function Location(interactions: Interactions, structure: Structure, unitA?: Unit, indexA?: Features.FeatureIndex, unitB?: Unit, indexB?: Features.FeatureIndex): Location {
        return DataLocation('interactions', { structure, interactions },
            { unitA: unitA as any, indexA: indexA as any, unitB: unitB as any, indexB: indexB as any });
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === 'interactions';
    }

    export function areLocationsEqual(locA: Location, locB: Location) {
        return (
            locA.data.structure === locB.data.structure &&
            locA.data.interactions === locB.data.interactions &&
            locA.element.indexA === locB.element.indexA &&
            locA.element.indexB === locB.element.indexB &&
            locA.element.unitA === locB.element.unitA &&
            locA.element.unitB === locB.element.unitB
        );
    }

    function _label(interactions: Interactions, element: Element): string {
        const { unitA, indexA, unitB, indexB } = element;
        const { contacts, unitsContacts } = interactions;
        if (unitA === unitB) {
            const contacts = unitsContacts.get(unitA.id);
            const idx = contacts.getDirectedEdgeIndex(indexA, indexB);
            return interactionTypeLabel(contacts.edgeProps.type[idx]);
        } else {
            const idx = contacts.getEdgeIndex(indexA, unitA, indexB, unitB);
            return interactionTypeLabel(contacts.edges[idx].props.type);
        }
    }

    export function locationLabel(location: Location): string {
        return _label(location.data.interactions, location.element);
    }

    export interface Loci extends DataLoci<StructureInteractions, Element> { }

    export function Loci(structure: Structure, interactions: Interactions, elements: ReadonlyArray<Element>): Loci {
        return DataLoci('interactions', { structure, interactions }, elements,
            (boundingSphere) => getBoundingSphere(interactions, elements, boundingSphere),
            () => getLabel(structure, interactions, elements));
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === 'interactions';
    }

    export function getBoundingSphere(interactions: Interactions, elements: ReadonlyArray<Element>, boundingSphere: Sphere3D) {
        const { unitsFeatures } = interactions;
        return CentroidHelper.fromPairProvider(elements.length, (i, pA, pB) => {
            const e = elements[i];
            Features.setPosition(pA, e.unitA, e.indexA, unitsFeatures.get(e.unitA.id));
            Features.setPosition(pB, e.unitB, e.indexB, unitsFeatures.get(e.unitB.id));
        }, boundingSphere);
    }

    export function getLabel(structure: Structure, interactions: Interactions, elements: ReadonlyArray<Element>) {
        const element = elements[0];
        if (element === undefined) return '';
        const { unitA, indexA, unitB, indexB } = element;
        const { unitsFeatures } = interactions;
        const { members: mA, offsets: oA } = unitsFeatures.get(unitA.id);
        const { members: mB, offsets: oB } = unitsFeatures.get(unitB.id);
        const options = { granularity: 'element' as LabelGranularity };
        if (oA[indexA + 1] - oA[indexA] > 1 || oB[indexB + 1] - oB[indexB] > 1) {
            options.granularity = 'residue';
        }
        return [
            _label(interactions, element),
            bondLabel(Bond.Location(structure, unitA, mA[oA[indexA]], structure, unitB, mB[oB[indexB]]), options)
        ].join('</br>');
    }
}

const FeatureProviders = [
    HydrogenDonorProvider, WeakHydrogenDonorProvider, HydrogenAcceptorProvider,
    NegativChargeProvider, PositiveChargeProvider, AromaticRingProvider,
    HalogenDonorProvider, HalogenAcceptorProvider,
    HydrophobicAtomProvider,
    MetalProvider, MetalBindingProvider,
];

const ContactProviders = {
    'ionic': IonicProvider,
    'pi-stacking': PiStackingProvider,
    'cation-pi': CationPiProvider,
    'halogen-bonds': HalogenBondsProvider,
    'hydrogen-bonds': HydrogenBondsProvider,
    'weak-hydrogen-bonds': WeakHydrogenBondsProvider,
    'hydrophobic': HydrophobicProvider,
    'metal-coordination': MetalCoordinationProvider,
};
type ContactProviders = typeof ContactProviders

function getProvidersParams(defaultOn: string[] = []) {
    const params: { [k in keyof ContactProviders]: PD.Mapped<PD.NamedParamUnion<{
        on: PD.Group<ContactProviders[k]['params']>
        off: PD.Group<{}>
    }>> } = Object.create(null);

    Object.keys(ContactProviders).forEach(k => {
        (params as any)[k] = PD.MappedStatic(defaultOn.includes(k) ? 'on' : 'off', {
            on: PD.Group(ContactProviders[k as keyof ContactProviders].params),
            off: PD.Group({})
        }, { cycle: true });
    });
    return params;
}
export const ContactProviderParams = getProvidersParams([
    // 'ionic',
    'cation-pi',
    'pi-stacking',
    'hydrogen-bonds',
    'halogen-bonds',
    // 'hydrophobic',
    'metal-coordination',
    // 'weak-hydrogen-bonds',
]);

export const InteractionsParams = {
    providers: PD.Group(ContactProviderParams, { isFlat: true }),
    contacts: PD.Group(ContactsParams, { label: 'Advanced Options' }),
};
export type InteractionsParams = typeof InteractionsParams
export type InteractionsProps = PD.Values<InteractionsParams>

export async function computeInteractions(ctx: CustomProperty.Context, structure: Structure, props: Partial<InteractionsProps>): Promise<Interactions> {
    const p = { ...PD.getDefaultValues(InteractionsParams), ...props };
    await ValenceModelProvider.attach(ctx, structure);

    const contactTesters: ContactTester[] = [];
    ObjectKeys(ContactProviders).forEach(k => {
        const { name, params } = p.providers[k];
        if (name === 'on') {
            contactTesters.push(ContactProviders[k].createTester(params as any));
        }
    });

    const requiredFeatures = new Set<FeatureType>();
    contactTesters.forEach(l => SetUtils.add(requiredFeatures, l.requiredFeatures));
    const featureProviders = FeatureProviders.filter(f => SetUtils.areIntersecting(requiredFeatures, f.types));

    const unitsFeatures = IntMap.Mutable<Features>();
    const unitsContacts = IntMap.Mutable<InteractionsIntraContacts>();

    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const group = structure.unitSymmetryGroups[i];
        if (ctx.runtime.shouldUpdate) {
            await ctx.runtime.update({ message: 'computing interactions', current: i, max: il });
        }
        const features = findUnitFeatures(structure, group.units[0], featureProviders);
        const intraUnitContacts = findIntraUnitContacts(structure, group.units[0], features, contactTesters, p.contacts);
        for (let j = 0, jl = group.units.length; j < jl; ++j) {
            const u = group.units[j];
            unitsFeatures.set(u.id, features);
            unitsContacts.set(u.id, intraUnitContacts);
        }
    }

    const contacts = findInterUnitContacts(structure, unitsFeatures, contactTesters, p.contacts);

    const interactions = { unitsFeatures, unitsContacts, contacts };
    refineInteractions(structure, interactions);
    return interactions;
}

function findUnitFeatures(structure: Structure, unit: Unit, featureProviders: Features.Provider[]) {
    const count = unit.elements.length;
    const featuresBuilder = FeaturesBuilder.create(count, count / 2);
    if (Unit.isAtomic(unit)) {
        for (const fp of featureProviders) {
            fp.add(structure, unit, featuresBuilder);
        }
    }
    return featuresBuilder.getFeatures(count);
}

function findIntraUnitContacts(structure: Structure, unit: Unit, features: Features, contactTesters: ReadonlyArray<ContactTester>, props: ContactsProps) {
    const builder = IntraContactsBuilder.create(features, unit.elements.length);
    if (Unit.isAtomic(unit)) {
        addUnitContacts(structure, unit, features, builder, contactTesters, props);
    }
    return builder.getContacts();
}

function findInterUnitContacts(structure: Structure, unitsFeatures: IntMap<Features>, contactTesters: ReadonlyArray<ContactTester>, props: ContactsProps) {
    const builder = InterContactsBuilder.create();

    Structure.eachUnitPair(structure, (unitA: Unit, unitB: Unit) => {
        const featuresA = unitsFeatures.get(unitA.id);
        const featuresB = unitsFeatures.get(unitB.id);
        addStructureContacts(structure, unitA as Unit.Atomic, featuresA, unitB as Unit.Atomic, featuresB, builder, contactTesters, props);
    }, {
        maxRadius: Math.max(...contactTesters.map(t => t.maxDistance)),
        validUnit: (unit: Unit) => Unit.isAtomic(unit),
        validUnitPair: (unitA: Unit, unitB: Unit) => Structure.validUnitPair(structure, unitA, unitB)
    });

    return builder.getContacts(unitsFeatures);
}