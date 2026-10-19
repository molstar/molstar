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
  for (const tester of testers) {
    _addUnitContacts(structure, unit, features, builder, tester, props);
  }
}

function _addUnitContacts(
  structure: Structure,
  unit: Unit.Atomic,
  features: Features,
  builder: IntraContactsBuilder,
  tester: ContactTester,
  props: ContactsProps,
) {
  const { x, y, z } = features;
  const { lookup3d, indices: subsetIndices } = features.subset(tester.requiredFeatures);

  const infoA = Features.Info(structure, unit, features);
  const infoB = { ...infoA };

  const distFactor = props.lineOfSightDistFactor;

  for (let t = 0, tl = OrderedSet.size(subsetIndices); t < tl; ++t) {
    const i = OrderedSet.getAt(subsetIndices, t);
    const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], tester.maxDistance);
    if (count === 0) continue;

    infoA.feature = i;

    for (let r = 0; r < count; ++r) {
      const j = OrderedSet.getAt(subsetIndices, indices[r]);
      if (j <= i) continue;

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

const _pairCandidates = PairContactCandidates();

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
  findStructureContactCandidates(structure, unitA, featuresA, unitB, featuresB, testers, _pairCandidates);
  addStructureContactCandidates(structure, unitA, featuresA, unitB, featuresB, _pairCandidates, builder, props);
}

/** The contacts of one unit pair before line of sight: feature pairs (i, j) and the first tester's type, in add order */
export interface PairContactCandidates {
  count: number;
  a: Int32Array;
  b: Int32Array;
  type: Int8Array;
}

export function PairContactCandidates(): PairContactCandidates {
  return { count: 0, a: new Int32Array(16), b: new Int32Array(16), type: new Int8Array(16) };
}

function pushCandidate(c: PairContactCandidates, i: number, j: number, type: InteractionType) {
  if (c.count === c.a.length) {
    const n = 2 * c.a.length;
    const a = new Int32Array(n);
    a.set(c.a);
    c.a = a;
    const b = new Int32Array(n);
    b.set(c.b);
    c.b = b;
    const t = new Int8Array(n);
    t.set(c.type);
    c.type = t;
  }
  c.a[c.count] = i;
  c.b[c.count] = j;
  c.type[c.count] = type;
  c.count++;
}

/**
 * The candidates of `addStructureContacts`: the same search, `validPair` and testers, but each feature pair with a type
 * is recorded instead of checked for line of sight and added. `addStructureContacts` adds a pair with the first
 * tester's type iff its line of sight is free (line of sight does not depend on the type, so a later tester never
 * adds it), so `addStructureContactCandidates` on these gives the same contacts in the same order.
 */
export function findStructureContactCandidates(
  structure: Structure,
  unitA: Unit.Atomic,
  featuresA: Features,
  unitB: Unit.Atomic,
  featuresB: Features,
  testers: ReadonlyArray<ContactTester>,
  out: PairContactCandidates,
) {
  out.count = 0;
  const { x: xA, y: yA, z: zA } = featuresA;
  const { lookup3d } = featuresB;

  const imageTransform = Mat4.mul(
    _imageTransform,
    unitB.conformation.operator.inverse,
    unitA.conformation.operator.matrix,
  );
  const isNotIdentity = !Mat4.isIdentity(imageTransform);
  const imageA = Vec3();

  const maxDistance = Math.max(...testers.map((t) => t.maxDistance));
  const { center, radius } = lookup3d.boundary.sphere;
  const testDistanceSq = (radius + maxDistance) * (radius + maxDistance);

  const infoA = Features.Info(structure, unitA, featuresA);
  const infoB = Features.Info(structure, unitB, featuresB);

  const candidates = getCandidates(featuresA, unitA, unitB, lookup3d.boundary.sphere, radius + maxDistance);
  for (let ci = 0, cl = candidates.length; ci < cl; ++ci) {
    const i = candidates[ci] as Features.FeatureIndex;
    Vec3.set(imageA, xA[i], yA[i], zA[i]);
    if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform);
    if (Vec3.squaredDistance(imageA, center) > testDistanceSq) continue;

    const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], maxDistance);
    if (count === 0) continue;

    infoA.feature = i;

    for (let r = 0; r < count; ++r) {
      const j = indices[r];
      infoB.feature = j;
      if (!validPair(structure, infoA, infoB)) continue;

      const distanceSq = squaredDistances[r];
      for (const tester of testers) {
        if (distanceSq < tester.maxDistance * tester.maxDistance) {
          const type = tester.getType(structure, infoA, infoB, distanceSq);
          if (type) {
            pushCandidate(out, i, j, type);
            break;
          }
        }
      }
    }
  }
}

/** Adds the candidates of a unit pair (from `findStructureContactCandidates`) whose line of sight is free */
export function addStructureContactCandidates(
  structure: Structure,
  unitA: Unit.Atomic,
  featuresA: Features,
  unitB: Unit.Atomic,
  featuresB: Features,
  candidates: PairContactCandidates,
  builder: InterContactsBuilder,
  props: ContactsProps,
) {
  const distFactor = props.lineOfSightDistFactor;
  const infoA = Features.Info(structure, unitA, featuresA);
  const infoB = Features.Info(structure, unitB, featuresB);
  const { a, b, type } = candidates;

  builder.startUnitPair(unitA, unitB);
  for (let k = 0, kl = candidates.count; k < kl; ++k) {
    infoA.feature = a[k] as Features.FeatureIndex;
    infoB.feature = b[k] as Features.FeatureIndex;
    if (checkLineOfSight(structure, infoA, infoB, distFactor)) {
      builder.add(infoA.feature, infoB.feature, type[k] as InteractionType);
    }
  }
  builder.finishUnitPair();
}
