/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Sphere3D } from '@molstar/core/math/geometry';
import { Result } from '@molstar/core/math/geometry/lookup3d/common';
import { type Structure, type StructureElement, Unit } from '@molstar/model/model/structure';
import { VdwRadius } from '@molstar/model/model/structure/model/properties/atomic';
import { Elements } from '@molstar/model/model/structure/model/properties/atomic/types';
import { StructureLookup3DResultContext } from '@molstar/model/model/structure/structure/util/lookup3d';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { altLoc, connectedTo, typeSymbol } from '../chemistry/util.js';
import type { FeatureType, InteractionType } from './common.js';
import type { InterContactsBuilder, IntraContactsBuilder } from './contacts-builder.js';
import { Features } from './features.js';
import type { ContactsProps } from './params.js';

const MAX_LINE_OF_SIGHT_DISTANCE = 3;

export interface ContactProvider<P extends PD.Params> {
  readonly name: string;
  readonly params: P;
  createTester(props: PD.Values<P>): ContactTester;
}

export interface ContactTester {
  readonly maxDistance: number;
  readonly requiredFeatures: ReadonlySet<FeatureType>;
  /**
   * The feature type pairs `getType` may accept: it returns `undefined` for any `(infoA, infoB)` whose types fail this
   * test. Without it, any two of `requiredFeatures`. The contact loops skip the other pairs before `validPair`, and
   * query a feature no farther than its type's partners need (perf-audit-2 06/P4, 06/P10).
   */
  readonly acceptsTypes?: (typeA: FeatureType, typeB: FeatureType) => boolean;
  getType: (
    structure: Structure,
    infoA: Features.Info,
    infoB: Features.Info,
    distanceSq: number,
  ) => InteractionType | undefined;
}

function validPair(structure: Structure, infoA: Features.Info, infoB: Features.Info): boolean {
  const indexA = infoA.members[infoA.offsets[infoA.feature]];
  const indexB = infoB.members[infoB.offsets[infoB.feature]];
  if (indexA === indexB && infoA.unit === infoB.unit) return false; // no self interaction

  const altA = altLoc(infoA.unit, indexA);
  const altB = altLoc(infoB.unit, indexB);
  if (altA && altB && altA !== altB) return false; // incompatible alternate location id
  if (
    infoA.unit === infoB.unit &&
    infoA.unit.model.atomicHierarchy.residueAtomSegments.count > 1 &&
    infoA.unit.residueIndex[infoA.unit.elements[indexA]] === infoB.unit.residueIndex[infoB.unit.elements[indexB]]
  )
    return false; // same residue (and more than one residue)

  // e.g. no hbond if donor and acceptor are bonded
  if (connectedTo(structure, infoA.unit, indexA, infoB.unit, indexB)) return false;

  return true;
}

//

function invalidAltLoc(
  unitA: Unit.Atomic,
  indexA: StructureElement.UnitIndex,
  unitB: Unit.Atomic,
  indexB: StructureElement.UnitIndex,
) {
  const altA = altLoc(unitA, indexA);
  const altB = altLoc(unitB, indexB);
  return altA && altB && altA !== altB;
}

function isMember(element: StructureElement.UnitIndex, info: Features.Info) {
  const { feature, offsets, members } = info;
  for (let i = offsets[feature], il = offsets[feature + 1]; i < il; ++i) {
    if (members[i] === element) return true;
  }
  return false;
}

const tmpVec = Vec3();
const tmpVecA = Vec3();
const tmpVecB = Vec3();

// need to use a separate context for structure.lookup3d.find because of nested queries
const lineOfSightLookupCtx = StructureLookup3DResultContext();

function checkLineOfSight(structure: Structure, infoA: Features.Info, infoB: Features.Info, distFactor: number) {
  const featureA = infoA.feature;
  const featureB = infoB.feature;
  const indexA = infoA.members[infoA.offsets[featureA]];
  const indexB = infoB.members[infoB.offsets[featureB]];

  Features.position(tmpVecA, infoA);
  Features.position(tmpVecB, infoB);
  Vec3.scale(tmpVec, Vec3.add(tmpVec, tmpVecA, tmpVecB), 0.5);

  const distMax = distFactor * MAX_LINE_OF_SIGHT_DISTANCE;

  const { count, indices, units, squaredDistances } = structure.lookup3d.find(
    tmpVec[0],
    tmpVec[1],
    tmpVec[2],
    distMax,
    lineOfSightLookupCtx,
  );
  if (count === 0) return true;

  for (let r = 0; r < count; ++r) {
    const i = indices[r];
    const unit = units[r];
    if (!Unit.isAtomic(unit)) continue;

    const element = typeSymbol(unit, i);
    // allow hydrogens
    if (element === Elements.H) continue;

    const vdw = VdwRadius(element);
    // check distance
    if (vdw * vdw * distFactor * distFactor <= squaredDistances[r]) continue;

    // allow different altlocs
    if (invalidAltLoc(unit, i, infoA.unit, indexA) || invalidAltLoc(unit, i, infoB.unit, indexB)) continue;

    // allow member atoms
    if ((infoA.unit === unit && isMember(i, infoA)) || (infoB.unit === unit && isMember(i, infoB))) continue;

    unit.conformation.position(unit.elements[i], tmpVec);
    // allow atoms at the center of functional groups
    if (Vec3.squaredDistance(tmpVec, tmpVecA) < 1 || Vec3.squaredDistance(tmpVec, tmpVecB) < 1) continue;

    return false;
  }

  return true;
}

/** Feature types are below 16 (`FeatureType`): a type set is a 16-bit mask, a type pair `typeA * 16 + typeB`. */
const TypeCount = 16;

/** The testers of one computation by the feature type pairs they accept (`ContactTester.acceptsTypes`) */
interface TesterTable {
  /** per tester: for each type, the mask of the types it accepts with it */
  readonly partners: readonly Uint16Array[];
  /** per type pair: the testers accepting it, in order */
  readonly pairTesters: readonly (readonly ContactTester[])[];
  /** per type pair: the longest `maxDistance` of its testers (0: none) */
  readonly pairDistance: Float64Array;
}

const testerTables = new WeakMap<ReadonlyArray<ContactTester>, TesterTable>();

function getTesterTable(testers: ReadonlyArray<ContactTester>) {
  let table = testerTables.get(testers);
  if (table) return table;
  const partners: Uint16Array[] = [];
  const pairTesters: ContactTester[][] = [];
  const pairDistance = new Float64Array(TypeCount * TypeCount);
  for (let p = 0; p < TypeCount * TypeCount; p++) pairTesters.push([]);
  for (const tester of testers) {
    const { requiredFeatures: required, acceptsTypes } = tester;
    const accepts = acceptsTypes ?? ((a: FeatureType, b: FeatureType) => required.has(a) && required.has(b));
    const masks = new Uint16Array(TypeCount);
    for (let a = 0; a < TypeCount; a++) {
      for (let b = 0; b < TypeCount; b++) {
        if (!accepts(a as FeatureType, b as FeatureType)) continue;
        masks[a] |= 1 << b;
        const p = a * TypeCount + b;
        pairTesters[p].push(tester);
        pairDistance[p] = Math.max(pairDistance[p], tester.maxDistance);
      }
    }
    partners.push(masks);
  }
  table = { partners, pairTesters, pairDistance };
  testerTables.set(testers, table);
  return table;
}

const featureTypeMasks = new WeakMap<Features, number>();

/** The mask of the types of `features` */
function featureTypes(features: Features) {
  let mask = featureTypeMasks.get(features);
  if (mask !== undefined) return mask;
  mask = 0;
  const { types, count } = features;
  for (let i = 0; i < count; i++) mask |= 1 << types[i];
  featureTypeMasks.set(features, mask);
  return mask;
}

/**
 * Add all intra-unit contacts, i.e. pairs of features
 */
export function addUnitContacts(
  structure: Structure,
  unit: Unit.Atomic,
  features: Features,
  builder: IntraContactsBuilder,
  testers: ReadonlyArray<ContactTester>,
  props: ContactsProps,
) {
  const { partners } = getTesterTable(testers);
  const present = featureTypes(features);
  for (let k = 0; k < testers.length; k++) {
    // a tester with no accepted pair of the unit's types adds nothing
    const masks = partners[k];
    let partnered = false;
    for (let t = 0; t < TypeCount; t++) {
      if (present & (1 << t) && masks[t] & present) {
        partnered = true;
        break;
      }
    }
    if (partnered) _addUnitContacts(structure, unit, features, builder, testers[k], masks, props);
  }
}

function _addUnitContacts(
  structure: Structure,
  unit: Unit.Atomic,
  features: Features,
  builder: IntraContactsBuilder,
  tester: ContactTester,
  partners: Uint16Array,
  props: ContactsProps,
) {
  const { x, y, z, types } = features;
  const present = featureTypes(features);
  const { lookup3d, indices: subsetIndices } = features.subset(tester.requiredFeatures);

  const infoA = Features.Info(structure, unit, features);
  const infoB = { ...infoA };

  const distFactor = props.lineOfSightDistFactor;

  for (let t = 0, tl = OrderedSet.size(subsetIndices); t < tl; ++t) {
    const i = OrderedSet.getAt(subsetIndices, t);
    const partnersA = partners[types[i]];
    if ((partnersA & present) === 0) continue;
    const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], tester.maxDistance);
    if (count === 0) continue;

    infoA.feature = i;

    for (let r = 0; r < count; ++r) {
      const j = OrderedSet.getAt(subsetIndices, indices[r]);
      if (j <= i || (partnersA & (1 << types[j])) === 0) continue;

      infoB.feature = j;
      if (!validPair(structure, infoA, infoB)) continue;

      const type = tester.getType(structure, infoA, infoB, squaredDistances[r]);
      if (type && checkLineOfSight(structure, infoA, infoB, distFactor)) {
        builder.add(i, j, type);
      }
    }
  }
}

const _imageTransform = Mat4();
/** per type of A: the longest distance of a tester accepting it with a type B has */
const _queryRadius = new Float64Array(TypeCount);
const _bToA = Mat4();
const _bCenterInA = Vec3();
// not the lookup's own result: symmetry copies share one Features object, so the inner query on B's lookup would
// overwrite it
const _candidateResult = Result.create<Features.FeatureIndex>();
let _candidates = new Int32Array(1024);

/** All features of `featuresA`, in order. */
function allCandidates(featuresA: Features) {
  const n = featuresA.count;
  if (_candidates.length < n) _candidates = new Int32Array(Math.max(n, 2 * _candidates.length));
  for (let i = 0; i < n; i++) _candidates[i] = i;
  return _candidates.subarray(0, n);
}

/**
 * The features of `unitA` (sorted) that may lie within `radius` of `sphereB` (in `unitB`'s frame): fetched from
 * `featuresA`'s grid around the sphere centre mapped into `unitA`'s frame, instead of transforming every feature of
 * `unitA` (06-model-props-script N1). The radius is padded so that the exact test in `unitB`'s frame still decides.
 * Every feature is a candidate, without a grid query (or building `featuresA`'s grid), when `unitA`'s boundary sphere
 * (which holds its features: atoms and centroids of atoms) lies inside that radius, mostly ions and waters next to a
 * chain, and for operators that do not keep distances (mmCIF allows them with a warning).
 */
function getCandidates(featuresA: Features, unitA: Unit, unitB: Unit, sphereB: Sphere3D, radius: number) {
  if (featuresA.count === 0) return _candidates.subarray(0, 0);
  Mat4.mul(_bToA, unitA.conformation.operator.inverse, unitB.conformation.operator.matrix);
  if (!Mat4.isRigid(_bToA, 1e-4)) return allCandidates(featuresA);
  Vec3.transformMat4(_bCenterInA, sphereB.center, _bToA);
  const sphereA = unitA.boundary.sphere;
  if (Vec3.distance(_bCenterInA, sphereA.center) + sphereA.radius <= radius) return allCandidates(featuresA);
  const { indices, count } = featuresA.lookup3d.find(
    _bCenterInA[0],
    _bCenterInA[1],
    _bCenterInA[2],
    radius * (1 + 1e-3) + 1e-3,
    _candidateResult,
  );
  if (_candidates.length < count) _candidates = new Int32Array(Math.max(count, 2 * _candidates.length));
  for (let i = 0; i < count; i++) _candidates[i] = indices[i];
  // ascending, as the scan was: keeps the contact order
  return _candidates.subarray(0, count).sort();
}

/**
 * Add all inter-unit contacts, i.e. pairs of features
 */
export function addStructureContacts(
  structure: Structure,
  unitA: Unit.Atomic,
  featuresA: Features,
  unitB: Unit.Atomic,
  featuresB: Features,
  builder: InterContactsBuilder,
  testers: ReadonlyArray<ContactTester>,
  props: ContactsProps,
) {
  const { x: xA, y: yA, z: zA, types: typesA } = featuresA;
  const { lookup3d, types: typesB } = featuresB;

  // query each feature of A only as far as a tester needs for its type and a type in B, skip types with no partner
  // there; then the hits only of accepted type pairs, before `validPair` (perf-audit-2 06/P4, 06/P10)
  const { pairTesters, pairDistance } = getTesterTable(testers);
  const presentA = featureTypes(featuresA);
  const presentB = featureTypes(featuresB);
  let maxDistance = 0;
  for (let a = 0; a < TypeCount; a++) {
    let r = 0;
    if (presentA & (1 << a)) {
      for (let b = 0; b < TypeCount; b++) {
        if (presentB & (1 << b)) r = Math.max(r, pairDistance[a * TypeCount + b]);
      }
    }
    _queryRadius[a] = r;
    maxDistance = Math.max(maxDistance, r);
  }
  if (maxDistance === 0) return;

  // the lookup queries need to happen in the "unitB space".
  // that means imageA = inverseOperB(operA(i))
  const imageTransform = Mat4.mul(
    _imageTransform,
    unitB.conformation.operator.inverse,
    unitA.conformation.operator.matrix,
  );
  const isNotIdentity = !Mat4.isIdentity(imageTransform);
  const imageA = Vec3();

  const { center, radius } = lookup3d.boundary.sphere;
  const testDistanceSq = (radius + maxDistance) * (radius + maxDistance);

  const distFactor = props.lineOfSightDistFactor;

  const infoA = Features.Info(structure, unitA, featuresA);
  const infoB = Features.Info(structure, unitB, featuresB);

  builder.startUnitPair(unitA, unitB);

  const candidates = getCandidates(featuresA, unitA, unitB, lookup3d.boundary.sphere, radius + maxDistance);
  for (let ci = 0, cl = candidates.length; ci < cl; ++ci) {
    const i = candidates[ci] as Features.FeatureIndex;
    const queryRadius = _queryRadius[typesA[i]];
    if (queryRadius === 0) continue;
    Vec3.set(imageA, xA[i], yA[i], zA[i]);
    if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform);
    if (Vec3.squaredDistance(imageA, center) > testDistanceSq) continue;

    // a smaller radius finds a subsequence of the hits, in the same order (the grid's cells and buckets)
    const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], queryRadius);
    if (count === 0) continue;

    infoA.feature = i;
    const pairA = typesA[i] * TypeCount;

    for (let r = 0; r < count; ++r) {
      const j = indices[r];
      const pair = pairA + typesB[j];
      const pairMax = pairDistance[pair];
      const distanceSq = squaredDistances[r];
      if (distanceSq >= pairMax * pairMax) continue;

      infoB.feature = j;
      if (!validPair(structure, infoA, infoB)) continue;

      for (const tester of pairTesters[pair]) {
        if (distanceSq < tester.maxDistance * tester.maxDistance) {
          const type = tester.getType(structure, infoA, infoB, distanceSq);
          if (type && checkLineOfSight(structure, infoA, infoB, distFactor)) {
            builder.add(i, j, type);
            break;
          }
        }
      }
    }
  }

  builder.finishUnitPair();
}
