/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, Mat4 } from '@molstar/core/math/linear-algebra';
import { MeshBuilder } from '../mesh-builder.js';
import { type Primitive, transformPrimitive } from '../../../primitive/primitive.js';
import { Cylinder, type CylinderProps, DefaultCylinderProps } from '../../../primitive/cylinder.js';
import { Prism } from '../../../primitive/prism.js';
import { polygon } from '../../../primitive/polygon.js';
import { hashFnv32a } from '@molstar/core/data/util';
import type { Ray3D } from '@molstar/core/math/geometry/primitives/ray3d';

/**
 * Cylinder primitives by their (rounded) props: entries of a hash with the props checked, as two props sets can share
 * a 32-bit hash; the last `MaxCachedCylinders` made only, as continuous radii (the clip gizmos while dragged) would
 * keep one per value (perf-audit-2 08/P9)
 */
const cylinderMap = new Map<number, { values: Int32Array; primitive: Primitive }[]>();
let cachedCylinderCount = 0;
export const MaxCachedCylinders = 256;
const up = Vec3.create(0, 1, 0);

const tmpCylinderDir = Vec3();
const tmpCylinderMatDir = Vec3();
const tmpCylinderCenter = Vec3();
const tmpCylinderMat = Mat4();
const tmpCylinderMatRot = Mat4();
const tmpCylinderScale = Vec3();
const tmpCylinderStart = Vec3();
const tmpUp = Vec3();

function setCylinderMat(m: Mat4, start: Vec3, dir: Vec3, length: number, matchDir: boolean) {
  Vec3.setMagnitude(tmpCylinderMatDir, dir, length / 2);
  Vec3.add(tmpCylinderCenter, start, tmpCylinderMatDir);
  // ensure the direction used to create the rotation is always pointing in the same
  // direction so the triangles of adjacent cylinder will line up
  if (matchDir) Vec3.matchDirection(tmpUp, up, tmpCylinderMatDir);
  else Vec3.copy(tmpUp, up);
  Vec3.set(tmpCylinderScale, 1, length, 1);
  Vec3.makeRotation(tmpCylinderMatRot, tmpUp, tmpCylinderMatDir);
  Mat4.scale(m, tmpCylinderMatRot, tmpCylinderScale);
  return Mat4.setTranslation(m, tmpCylinderCenter);
}

const tmpPropValues = new Int32Array(9);
function getCylinderPropsKey(props: CylinderProps) {
  tmpPropValues[0] = Math.round(1000 * (props.radiusTop ?? DefaultCylinderProps.radiusTop));
  tmpPropValues[1] = Math.round(1000 * (props.radiusBottom ?? DefaultCylinderProps.radiusBottom));
  tmpPropValues[2] = Math.round(1000 * (props.height ?? DefaultCylinderProps.height));
  tmpPropValues[3] = props.radialSegments ?? DefaultCylinderProps.radialSegments;
  tmpPropValues[4] = props.heightSegments ?? DefaultCylinderProps.heightSegments;
  tmpPropValues[5] = (props.topCap ?? DefaultCylinderProps.topCap) ? 1 : 0;
  tmpPropValues[6] = (props.bottomCap ?? DefaultCylinderProps.bottomCap) ? 1 : 0;
  tmpPropValues[7] = Math.round(1000 * (props.thetaStart ?? DefaultCylinderProps.thetaStart));
  tmpPropValues[8] = Math.round(1000 * (props.thetaLength ?? DefaultCylinderProps.thetaLength));
  return hashFnv32a(tmpPropValues);
}

function sameValues(a: Int32Array, b: Int32Array) {
  for (let i = 0; i < 9; i++) if (a[i] !== b[i]) return false;
  return true;
}

function getCylinder(props: CylinderProps) {
  const key = getCylinderPropsKey(props);
  const entries = cylinderMap.get(key);
  if (entries) {
    for (const e of entries) {
      if (sameValues(e.values, tmpPropValues)) return e.primitive;
    }
  }
  let cylinder: Primitive;
  if (props.radialSegments && props.radialSegments <= 4) {
    const sideCount = Math.max(3, props.radialSegments);
    const prism = Prism(polygon(sideCount, true, props.radiusTop), props);
    cylinder = transformPrimitive(prism, Mat4.rotX90);
  } else {
    cylinder = Cylinder(props);
  }
  if (cachedCylinderCount >= MaxCachedCylinders) {
    // the hash first set (insertion order) goes; the cache is a reuse of shapes, not of a given one
    const first = cylinderMap.keys().next().value!;
    cachedCylinderCount -= cylinderMap.get(first)!.length;
    cylinderMap.delete(first);
  }
  const entry = { values: tmpPropValues.slice(), primitive: cylinder };
  if (entries && cylinderMap.has(key)) entries.push(entry);
  else cylinderMap.set(key, [entry]);
  cachedCylinderCount++;
  return cylinder;
}

export type BasicCylinderProps = Omit<CylinderProps, 'height'>;

export function addSimpleCylinder(state: MeshBuilder.State, start: Vec3, end: Vec3, props: BasicCylinderProps) {
  const d = Vec3.distance(start, end);
  Vec3.sub(tmpCylinderDir, end, start);
  setCylinderMat(tmpCylinderMat, start, tmpCylinderDir, d, false);
  MeshBuilder.addPrimitive(state, tmpCylinderMat, getCylinder(props));
}

export function addCylinderFromRay3D(state: MeshBuilder.State, ray: Ray3D, length: number, props: BasicCylinderProps) {
  setCylinderMat(tmpCylinderMat, ray.origin, ray.direction, length, false);
  MeshBuilder.addPrimitive(state, tmpCylinderMat, getCylinder(props));
}

export function addCylinder(
  state: MeshBuilder.State,
  start: Vec3,
  end: Vec3,
  lengthScale: number,
  props: BasicCylinderProps,
) {
  const d = Vec3.distance(start, end) * lengthScale;
  Vec3.sub(tmpCylinderDir, end, start);
  setCylinderMat(tmpCylinderMat, start, tmpCylinderDir, d, true);
  MeshBuilder.addPrimitive(state, tmpCylinderMat, getCylinder(props));
}

export function addDoubleCylinder(
  state: MeshBuilder.State,
  start: Vec3,
  end: Vec3,
  lengthScale: number,
  shift: Vec3,
  props: BasicCylinderProps,
) {
  const d = Vec3.distance(start, end) * lengthScale;
  const cylinder = getCylinder(props);
  Vec3.sub(tmpCylinderDir, end, start);
  // positivly shifted cylinder
  Vec3.add(tmpCylinderStart, start, shift);
  setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d, true);
  MeshBuilder.addPrimitive(state, tmpCylinderMat, cylinder);
  // negativly shifted cylinder
  Vec3.sub(tmpCylinderStart, start, shift);
  setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d, true);
  MeshBuilder.addPrimitive(state, tmpCylinderMat, cylinder);
}

export function addFixedCountDashedCylinder(
  state: MeshBuilder.State,
  start: Vec3,
  end: Vec3,
  lengthScale: number,
  segmentCount: number,
  stubCap: boolean,
  props: BasicCylinderProps,
) {
  const d = Vec3.distance(start, end) * lengthScale;
  const isOdd = segmentCount % 2 !== 0;
  const s = Math.floor((segmentCount + 1) / 2);
  let step = d / (segmentCount + 0.5);

  let cylinder = getCylinder(props);
  Vec3.setMagnitude(tmpCylinderDir, Vec3.sub(tmpCylinderDir, end, start), step);
  Vec3.copy(tmpCylinderStart, start);
  for (let j = 0; j < s; ++j) {
    Vec3.add(tmpCylinderStart, tmpCylinderStart, tmpCylinderDir);
    if (isOdd && j === s - 1) {
      if (!stubCap && props.topCap) {
        // a copy: the caller passes the same props again (MirroredAromatic's second row, perf-audit-2 08/B7)
        cylinder = getCylinder({ ...props, topCap: false });
      }
      step /= 2;
    }
    setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, step, false);
    MeshBuilder.addPrimitive(state, tmpCylinderMat, cylinder);

    Vec3.add(tmpCylinderStart, tmpCylinderStart, tmpCylinderDir);
  }
}
