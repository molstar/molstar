/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Structure, Unit, Bond, StructureElement } from '@molstar/model/model/structure';
import { Features, FeaturesBuilder } from './features.js';
import { ValenceModelProvider } from '../valence-model.js';
import {
  type InteractionsIntraContacts,
  type InteractionsInterContacts,
  type FeatureType,
  type InteractionType,
  type InteractionFlag,
  interactionTypeLabel,
} from './common.js';
import { IntraContactsBuilder, InterContactsBuilder } from './contacts-builder.js';
import { IntMap, OrderedSet } from '@molstar/core/data/int';
import {
  addStructureContactCandidates,
  addStructureContacts,
  addUnitContacts,
  type ContactTester,
  findStructureContactCandidates,
  PairContactCandidates,
} from './contacts.js';
import { BridgeProviderParams, InteractionsParams, type InteractionsProps, type ContactsProps } from './params.js';
import { HalogenDonorProvider, HalogenAcceptorProvider, HalogenBondsProvider } from './halogen-bonds.js';
import {
  HydrogenDonorProvider,
  WeakHydrogenDonorProvider,
  HydrogenAcceptorProvider,
  HydrogenBondsProvider,
  WeakHydrogenBondsProvider,
} from './hydrogen-bonds.js';
import { WaterBridgesProvider } from './water-bridges.js';
import {
  NegativChargeProvider,
  PositiveChargeProvider,
  AromaticRingProvider,
  IonicProvider,
  PiStackingProvider,
  CationPiProvider,
} from './charged.js';
import { HydrophobicAtomProvider, HydrophobicProvider } from './hydrophobic.js';
import { SetUtils } from '@molstar/core/util/set';
import { MetalCoordinationProvider, MetalProvider, MetalBindingProvider } from './metal.js';
import { refineInteractions } from './refine.js';
import type { CustomProperty } from '../../common/custom-property.js';
import { DataLocation } from '@molstar/model/model/location';
import { CentroidHelper } from '@molstar/core/math/geometry/centroid-helper';
import type { Sphere3D } from '@molstar/core/math/geometry';
import { DataLoci } from '@molstar/model/model/loci';
import { bondLabel, bundleLabel, type LabelGranularity } from '@molstar/model/model/label';
import { ObjectKeys } from '@molstar/core/util/type-helpers';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { arrayEqual } from '@molstar/core/util/array';
import { BondType } from '@molstar/model/model/structure/model/types';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';

export { Interactions, Bridges };
export type { BridgeContact, BridgeContacts };

interface BridgeContact {
  readonly unitA: number;
  readonly indexA: Features.FeatureIndex;
  readonly unitB: number;
  readonly indexB: Features.FeatureIndex;
  /** mediator unit id */
  readonly unitM: number;
  /** mediator feature facing endpoint A */
  readonly indexMA: Features.FeatureIndex;
  /** mediator feature facing endpoint B */
  readonly indexMB: Features.FeatureIndex;
  props: { type: InteractionType; flag: InteractionFlag };
}
type BridgeContacts = ReadonlyArray<BridgeContact>;

interface Interactions {
  /** Features of each unit */
  unitsFeatures: IntMap<Features>;
  /** Interactions of each unit */
  unitsContacts: IntMap<InteractionsIntraContacts>;
  /** Interactions between units */
  contacts: InteractionsInterContacts;
  /** Bridge-mediated interactions covering the whole structure */
  bridges: BridgeContacts;
}

namespace Interactions {
  type StructureInteractions = { readonly structure: Structure; readonly interactions: Interactions };

  export interface Element {
    unitA: Unit;
    /** Index into features of unitA */
    indexA: Features.FeatureIndex;
    unitB: Unit;
    /** Index into features of unitB */
    indexB: Features.FeatureIndex;
  }

  export interface Location extends DataLocation<StructureInteractions, Element> {}

  export function Location(
    interactions: Interactions,
    structure: Structure,
    unitA?: Unit,
    indexA?: Features.FeatureIndex,
    unitB?: Unit,
    indexB?: Features.FeatureIndex,
  ): Location {
    return DataLocation(
      'interactions',
      { structure, interactions },
      { unitA: unitA as any, indexA: indexA as any, unitB: unitB as any, indexB: indexB as any },
    );
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
      const idx = contacts.getEdgeIndex(indexA, unitA.id, indexB, unitB.id);
      return interactionTypeLabel(contacts.edges[idx].props.type);
    }
  }

  export function locationLabel(location: Location): string {
    return _label(location.data.interactions, location.element);
  }

  export interface Loci extends DataLoci<StructureInteractions, Element> {}

  export function Loci(structure: Structure, interactions: Interactions, elements: ReadonlyArray<Element>): Loci {
    return DataLoci(
      'interactions',
      { structure, interactions },
      elements,
      (boundingSphere) => getBoundingSphere(interactions, elements, boundingSphere),
      () => getLabel(structure, interactions, elements),
    );
  }

  export function isLoci(x: any): x is Loci {
    return !!x && x.kind === 'data-loci' && x.tag === 'interactions';
  }

  export function getBoundingSphere(
    interactions: Interactions,
    elements: ReadonlyArray<Element>,
    boundingSphere: Sphere3D,
  ) {
    const { unitsFeatures } = interactions;
    return CentroidHelper.fromPairProvider(
      elements.length,
      (i, pA, pB) => {
        const e = elements[i];
        Features.setPosition(pA, e.unitA, e.indexA, unitsFeatures.get(e.unitA.id));
        Features.setPosition(pB, e.unitB, e.indexB, unitsFeatures.get(e.unitB.id));
      },
      boundingSphere,
    );
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
      bondLabel(Bond.Location(structure, unitA, mA[oA[indexA]], structure, unitB, mB[oB[indexB]]), options),
    ].join('</br>');
  }
}

namespace Bridges {
  export interface Data {
    readonly structure: Structure;
    readonly bridges: BridgeContacts;
    readonly unitsFeatures: IntMap<Features>;
  }

  export interface Element {
    bridgeIndex: number;
  }

  export interface Location extends DataLocation<Data, Element> {}

  export function Location(data: Data, bridgeIndex = 0): Location {
    return DataLocation('bridges', data, { bridgeIndex });
  }

  export function isLocation(x: any): x is Location {
    return !!x && x.kind === 'data-location' && x.tag === 'bridges';
  }

  export interface Loci extends DataLoci<Data, Element> {}

  export function Loci(data: Data, elements: ReadonlyArray<Element>): Loci {
    return DataLoci(
      'bridges',
      data,
      elements,
      (bs) => getBoundingSphere(data, elements, bs),
      () => getLabel(data, elements),
    );
  }

  export function isLoci(x: any): x is Loci {
    return !!x && x.kind === 'data-loci' && x.tag === 'bridges';
  }

  function getLabel(data: Data, elements: ReadonlyArray<Element>): string {
    const e = elements[0];
    if (e === undefined) return '';

    const { structure, bridges, unitsFeatures } = data;
    const bridge = bridges[e.bridgeIndex];

    const uA = structure.unitMap.get(bridge.unitA) as Unit.Atomic;
    const fA = unitsFeatures.get(bridge.unitA);
    const uM = structure.unitMap.get(bridge.unitM) as Unit.Atomic;
    const fM = unitsFeatures.get(bridge.unitM);
    const uB = structure.unitMap.get(bridge.unitB) as Unit.Atomic;
    const fB = unitsFeatures.get(bridge.unitB);

    const options = { granularity: 'element' as LabelGranularity };
    if (
      fA.offsets[bridge.indexA + 1] - fA.offsets[bridge.indexA] > 1 ||
      fB.offsets[bridge.indexB + 1] - fB.offsets[bridge.indexB] > 1
    ) {
      options.granularity = 'residue';
    }

    return [
      interactionTypeLabel(bridge.props.type),
      bundleLabel(
        {
          loci: [
            StructureElement.Loci(structure, [
              {
                unit: uA,
                indices: OrderedSet.ofSingleton(fA.members[fA.offsets[bridge.indexA]] as StructureElement.UnitIndex),
              },
            ]),
            StructureElement.Loci(structure, [
              {
                unit: uM,
                indices: OrderedSet.ofSingleton(fM.members[fM.offsets[bridge.indexMA]] as StructureElement.UnitIndex),
              },
            ]),
            StructureElement.Loci(structure, [
              {
                unit: uB,
                indices: OrderedSet.ofSingleton(fB.members[fB.offsets[bridge.indexB]] as StructureElement.UnitIndex),
              },
            ]),
          ],
        },
        options,
      ),
    ].join('</br>');
  }

  function getBoundingSphere(data: Data, elements: ReadonlyArray<Element>, boundingSphere: Sphere3D) {
    return CentroidHelper.fromPairProvider(
      elements.length * 2,
      (i, pA, pB) => {
        const bridge = data.bridges[elements[i >> 1].bridgeIndex];

        const uA = data.structure.unitMap.get(bridge.unitA) as Unit.Atomic;
        const fA = data.unitsFeatures.get(bridge.unitA);
        const uM = data.structure.unitMap.get(bridge.unitM) as Unit.Atomic;
        const fM = data.unitsFeatures.get(bridge.unitM);
        const uB = data.structure.unitMap.get(bridge.unitB) as Unit.Atomic;
        const fB = data.unitsFeatures.get(bridge.unitB);

        const aIdx = fA.members[fA.offsets[bridge.indexA]];
        const mIdx = fM.members[fM.offsets[bridge.indexMA]];
        const bIdx = fB.members[fB.offsets[bridge.indexB]];

        if ((i & 1) === 0) {
          uA.conformation.position(uA.elements[aIdx], pA);
          uM.conformation.position(uM.elements[mIdx], pB);
        } else {
          uM.conformation.position(uM.elements[mIdx], pA);
          uB.conformation.position(uB.elements[bIdx], pB);
        }
      },
      boundingSphere,
    );
  }
}

const FeatureProviders = [
  HydrogenDonorProvider,
  WeakHydrogenDonorProvider,
  HydrogenAcceptorProvider,
  NegativChargeProvider,
  PositiveChargeProvider,
  AromaticRingProvider,
  HalogenDonorProvider,
  HalogenAcceptorProvider,
  HydrophobicAtomProvider,
  MetalProvider,
  MetalBindingProvider,
];

const ContactProviders = {
  ionic: IonicProvider,
  'pi-stacking': PiStackingProvider,
  'cation-pi': CationPiProvider,
  'halogen-bonds': HalogenBondsProvider,
  'hydrogen-bonds': HydrogenBondsProvider,
  'weak-hydrogen-bonds': WeakHydrogenBondsProvider,
  hydrophobic: HydrophobicProvider,
  'metal-coordination': MetalCoordinationProvider,
};
type ContactProviders = typeof ContactProviders;

const BridgeProviders = {
  'water-bridges': WaterBridgesProvider,
};
type BridgeProviders = typeof BridgeProviders;

export interface ComputeInterctionsOptions {
  skipIntraContacts?: boolean;
  unitPairTest?: (unitA: Unit, unitB: Unit) => boolean;
  /** Filled with the unit pairs of the inter-unit search */
  stats?: InterUnitContactsStats;
}

export interface InterUnitContactsStats {
  /** close unit pairs */
  pairs: number;
  /** pairs searched as the first of their symmetry-equivalent class ('per-pair-class') */
  classes: number;
  /** pairs that reused their class's search ('per-pair-class'); the rest were searched in full */
  replayed: number;
}

export async function computeInteractions(
  ctx: CustomProperty.Context,
  structure: Structure,
  props: Partial<InteractionsProps>,
  options?: ComputeInterctionsOptions,
): Promise<Interactions> {
  const p = { ...PD.getDefaultValues(InteractionsParams), ...props };
  // features do not depend on how unit pairs are searched
  const { symmetry: _, ...featureProps } = p;
  const cacheKey = JSON.stringify(featureProps);

  await ValenceModelProvider.attach(ctx, structure);

  const contactTesters: ContactTester[] = [];
  ObjectKeys(ContactProviders).forEach((k) => {
    const { name, params } = p.providers[k];
    if (name === 'on') {
      contactTesters.push(ContactProviders[k].createTester(params as any));
    }
  });

  const requiredFeatures = new Set<FeatureType>();
  contactTesters.forEach((l) => SetUtils.add(requiredFeatures, l.requiredFeatures));
  ObjectKeys(BridgeProviders).forEach((k) => {
    if (p.bridges[k].name === 'on') SetUtils.add(requiredFeatures, BridgeProviders[k].requiredFeatures);
  });
  const featureProviders = FeatureProviders.filter((f) => SetUtils.areIntersecting(requiredFeatures, f.types));

  const unitsFeatures = IntMap.Mutable<Features>();
  const unitsContacts = IntMap.Mutable<InteractionsIntraContacts>();

  for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
    const group = structure.unitSymmetryGroups[i];
    if (ctx.runtime.shouldUpdate) {
      await ctx.runtime.update({ message: 'computing interactions', current: i, max: il });
    }
    const features = findUnitFeatures(structure, group.units[0], featureProviders, cacheKey);
    for (let j = 0, jl = group.units.length; j < jl; ++j) {
      const u = group.units[j];
      unitsFeatures.set(u.id, features);
    }

    if (options?.skipIntraContacts) continue;

    const intraUnitContacts = findIntraUnitContacts(structure, group.units[0], features, contactTesters, p.contacts);
    for (let j = 0, jl = group.units.length; j < jl; ++j) {
      const u = group.units[j];
      unitsContacts.set(u.id, intraUnitContacts);
    }
  }

  const contacts =
    p.symmetry === 'per-pair-class'
      ? findInterUnitContactsPerPairClass(structure, unitsFeatures, contactTesters, p.contacts, options)
      : findInterUnitContacts(structure, unitsFeatures, contactTesters, p.contacts, options);
  const bridges = findBridges(structure, unitsFeatures, p.bridges);
  const interactions = { unitsFeatures, unitsContacts, contacts, bridges };

  refineInteractions(structure, interactions);
  return interactions;
}

function findUnitFeatures(structure: Structure, unit: Unit, featureProviders: Features.Provider[], cacheKey: string) {
  const key = `features-${cacheKey}`;
  if (unit.transientCache.has(key)) {
    return unit.transientCache.get(key) as Features;
  }

  const count = unit.elements.length;
  const featuresBuilder = FeaturesBuilder.create(count, count / 2);
  if (Unit.isAtomic(unit)) {
    for (const fp of featureProviders) {
      fp.add(structure, unit, featuresBuilder);
    }
  }
  const features = featuresBuilder.getFeatures(count);
  unit.transientCache.set(key, features);
  return features;
}

function findIntraUnitContacts(
  structure: Structure,
  unit: Unit,
  features: Features,
  contactTesters: ReadonlyArray<ContactTester>,
  props: ContactsProps,
) {
  const builder = IntraContactsBuilder.create(features, unit.elements.length);
  if (Unit.isAtomic(unit)) {
    addUnitContacts(structure, unit, features, builder, contactTesters, props);
  }
  return builder.getContacts();
}

function findBridges(
  structure: Structure,
  unitsFeatures: IntMap<Features>,
  props: PD.Values<typeof BridgeProviderParams>,
): BridgeContacts {
  const bridges: BridgeContact[] = [];

  ObjectKeys(BridgeProviders).forEach((k) => {
    const { name, params } = props[k];
    if (name === 'on') {
      for (const b of BridgeProviders[k].find(structure, unitsFeatures, params as any)) bridges.push(b);
    }
  });

  return bridges;
}

function findInterUnitContacts(
  structure: Structure,
  unitsFeatures: IntMap<Features>,
  contactTesters: ReadonlyArray<ContactTester>,
  props: ContactsProps,
  options?: ComputeInterctionsOptions,
) {
  const builder = InterContactsBuilder.create();

  const { stats } = options ?? {};
  Structure.eachUnitPair(
    structure,
    (unitA: Unit, unitB: Unit) => {
      if (stats) stats.pairs++;
      const featuresA = unitsFeatures.get(unitA.id);
      const featuresB = unitsFeatures.get(unitB.id);
      addStructureContacts(
        structure,
        unitA as Unit.Atomic,
        featuresA,
        unitB as Unit.Atomic,
        featuresB,
        builder,
        contactTesters,
        props,
      );
    },
    unitPairOptions(structure, contactTesters, options),
  );

  return builder.getContacts(unitsFeatures);
}

function unitPairOptions(
  structure: Structure,
  contactTesters: ReadonlyArray<ContactTester>,
  options?: ComputeInterctionsOptions,
): Structure.EachUnitPairProps {
  return {
    maxRadius: Math.max(...contactTesters.map((t) => t.maxDistance)),
    validUnit: (unit: Unit) => Unit.isAtomic(unit),
    validUnitPair: (unitA: Unit, unitB: Unit) => {
      if (options?.unitPairTest) return options.unitPairTest(unitA, unitB);
      return Structure.validUnitPair(structure, unitA, unitB);
    },
  };
}

/**
 * Units whose inter-unit contacts are a function of their symmetry group and operator:
 * - a rigid operator (the testers measure world distances and angles);
 * - no covalent bond to another unit (`validPair` and the testers' bond walks read covalent inter-unit bonds);
 * - intra-unit bonds not scoped by operator (`IndexPairBonds` with operators can bond units of one group differently);
 * - the same ideal atom geometries as the group's first unit that has none of the above (the valence model is per
 *   unit and comes from the root's topology, e.g. bonds a child structure does not have; the H-bond testers read it).
 */
function getPairClassUnits(structure: Structure) {
  const eligible = new Set<number>();
  const bonded = new Set<number>();
  for (const e of structure.interUnitBonds.edges) {
    if (BondType.isCovalent(e.props.flag)) {
      bonded.add(e.unitA);
      bonded.add(e.unitB);
    }
  }
  const valenceModel = ValenceModelProvider.get(structure).value!;
  for (const group of structure.unitSymmetryGroups) {
    const u0 = group.units[0];
    if (!Unit.isAtomic(u0) || IndexPairBonds.Provider.get(u0.model)?.hasOperators) continue;
    let reference: ArrayLike<number> | undefined;
    for (const u of group.units) {
      if (bonded.has(u.id) || !Mat4.isRigid(u.conformation.operator.matrix, PairClassRigidEps)) continue;
      const idealGeometry = valenceModel.get(u.id)?.idealGeometry;
      if (!idealGeometry) continue;
      reference ??= idealGeometry;
      if (arrayEqual(idealGeometry, reference)) eligible.add(u.id);
    }
  }
  return eligible;
}

interface PairClass {
  /** unitB.operator⁻¹ · unitA.operator */
  readonly transform: Mat4;
  readonly candidates: PairContactCandidates;
}

const _pairTransform = Mat4();
/**
 * Entries of the relative operator rounded for the class key; classes are then told apart by `Mat4.areEqual`. Two
 * equal operators rounded to different keys only make two classes (each searched once).
 */
const PairClassKeyScale = 1e3;
const PairClassTransformEps = 1e-7;
/** Operators printed with 6 decimals are rigid to about 1e-6 */
const PairClassRigidEps = 1e-5;

function pairClassKey(groupA: number, groupB: number, t: Mat4) {
  let key = `${groupA} ${groupB}`;
  for (let i = 0; i < 16; i++) key += ` ${Math.round(t[i] * PairClassKeyScale)}`;
  return key;
}

/**
 * Inter-unit contacts searched once per class of symmetry-equivalent unit pairs: the same two symmetry groups (so the
 * same features) and the same relative operator (so the same distances and angles). The first pair of a class runs
 * the search (`findStructureContactCandidates`) and the others reuse its candidates; line of sight, which sees third
 * units, is checked per pair. Pairs with a unit outside `getPairClassUnits` are searched in full. Edges are added in
 * the order of `findInterUnitContacts`.
 */
function findInterUnitContactsPerPairClass(
  structure: Structure,
  unitsFeatures: IntMap<Features>,
  contactTesters: ReadonlyArray<ContactTester>,
  props: ContactsProps,
  options?: ComputeInterctionsOptions,
) {
  const builder = InterContactsBuilder.create();
  const { stats } = options ?? {};
  const eligible = getPairClassUnits(structure);
  const groupIndex = structure.unitSymmetryGroupsIndexMap;
  const classes = new Map<string, PairClass[]>();
  const fullCandidates = PairContactCandidates();

  Structure.eachUnitPair(
    structure,
    (unitA: Unit, unitB: Unit) => {
      if (stats) stats.pairs++;
      const uA = unitA as Unit.Atomic;
      const uB = unitB as Unit.Atomic;
      const featuresA = unitsFeatures.get(unitA.id);
      const featuresB = unitsFeatures.get(unitB.id);

      if (!eligible.has(unitA.id) || !eligible.has(unitB.id)) {
        findStructureContactCandidates(structure, uA, featuresA, uB, featuresB, contactTesters, fullCandidates);
        addStructureContactCandidates(structure, uA, featuresA, uB, featuresB, fullCandidates, builder, props);
        return;
      }

      const t = Mat4.mul(_pairTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix);
      const key = pairClassKey(groupIndex.get(unitA.id), groupIndex.get(unitB.id), t);
      let bucket = classes.get(key);
      if (!bucket) classes.set(key, (bucket = []));
      let pairClass: PairClass | undefined;
      for (const c of bucket) {
        if (Mat4.areEqual(c.transform, t, PairClassTransformEps)) {
          pairClass = c;
          break;
        }
      }

      if (pairClass) {
        if (stats) stats.replayed++;
      } else {
        if (stats) stats.classes++;
        pairClass = { transform: Mat4.copy(Mat4(), t), candidates: PairContactCandidates() };
        bucket.push(pairClass);
        findStructureContactCandidates(structure, uA, featuresA, uB, featuresB, contactTesters, pairClass.candidates);
      }
      addStructureContactCandidates(structure, uA, featuresA, uB, featuresB, pairClass.candidates, builder, props);
    },
    unitPairOptions(structure, contactTesters, options),
  );

  return builder.getContacts(unitsFeatures);
}
