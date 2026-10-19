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

/** Whether the linear part of `m` is orthonormal within `eps` (a rigid transform keeps distances). */
function isRigid(m: Mat4, eps: number) {
  for (let i = 0; i < 3; i++) {
    for (let j = i; j < 3; j++) {
      const d = m[4 * i] * m[4 * j] + m[4 * i + 1] * m[4 * j + 1] + m[4 * i + 2] * m[4 * j + 2];
      if (Math.abs(d - (i === j ? 1 : 0)) > eps) return false;
    }
  }
  return true;
}

/**
 * The features of `unitA` (sorted) that may lie within `radius` of `sphereB` (in `unitB`'s frame): fetched from
 * `featuresA`'s grid around the sphere centre mapped into `unitA`'s frame, instead of transforming every feature of
 * `unitA` (06-model-props-script N1). The radius is padded so that the exact test in `unitB`'s frame still decides;
 * operators that do not keep distances (mmCIF allows them with a warning) get every feature, as before.
 */
function getCandidates(featuresA: Features, unitA: Unit, unitB: Unit, sphereB: Sphere3D, radius: number) {
  const n = featuresA.count;
  if (n === 0) return _candidates.subarray(0, 0);
  Mat4.mul(_bToA, unitA.conformation.operator.inverse, unitB.conformation.operator.matrix);
  if (!isRigid(_bToA, 1e-4)) {
    if (_candidates.length < n) _candidates = new Int32Array(Math.max(n, 2 * _candidates.length));
    for (let i = 0; i < n; i++) _candidates[i] = i;
    return _candidates.subarray(0, n);
  }
  Vec3.transformMat4(_bCenterInA, sphereB.center, _bToA);
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
  const { x: xA, y: yA, z: zA } = featuresA;
  const { lookup3d } = featuresB;

  // the lookup queries need to happen in the "unitB space".
  // that means imageA = inverseOperB(operA(i))
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

  const distFactor = props.lineOfSightDistFactor;

  const infoA = Features.Info(structure, unitA, featuresA);
  const infoB = Features.Info(structure, unitB, featuresB);

  builder.startUnitPair(unitA, unitB);

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
