import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { RuntimeContext } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ColorNames } from '@molstar/core/util/color/names';
import { ShapeRepresentation } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { Shape } from '@molstar/model/model/shape';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { structureElementLociLabelMany } from '@molstar/model/model/label';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { MarkerActions } from '@molstar/core/util/marker-action';
import { Plane } from '@molstar/graphics/geo/primitive/plane';
import { Circle } from '@molstar/graphics/geo/primitive/circle';
import { StructureElement } from '@molstar/model/model/structure';
import { Axes3D } from '@molstar/core/math/geometry';

export interface PlaneData {
  locis: StructureElement.Loci[];
}

const _PlaneParams = {
  ...Mesh.Params,
  color: PD.Color(ColorNames.orange),
  scaleFactor: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }),
  variant: PD.Select('rectangle', PD.arrayToOptions(['rectangle', 'circle'] as const), { label: 'Shape' }),
  radialSegments: PD.Numeric(64, { min: 12, max: 256, step: 1 }, { hideIf: (p) => p.variant !== 'circle' }),
};
type _PlaneParams = typeof _PlaneParams;

const PlaneVisuals = {
  plane: (ctx: RepresentationContext, getParams: RepresentationParamsGetter<PlaneData, _PlaneParams>) =>
    ShapeRepresentation(getPlaneShape, Mesh.Utils),
};

export const PlaneParams = {
  ..._PlaneParams,
  visuals: PD.MultiSelect(['plane'], PD.objectToOptions(PlaneVisuals)),
};
export type PlaneParams = typeof PlaneParams;
export type PlaneProps = PD.Values<PlaneParams>;

//

function getPlaneName(locis: StructureElement.Loci[]) {
  const label = structureElementLociLabelMany(locis, { countsOnly: true });
  return `Best Fit Plane of ${label}`;
}

const tmpMat = Mat4();
const tmpV = Vec3();
function buildPlaneMesh(data: PlaneData, props: PlaneProps, mesh?: Mesh): Mesh {
  const state = MeshBuilder.createState(256, 128, mesh);
  const principalAxes = StructureElement.Loci.getPrincipalAxesMany(data.locis);
  const axes = principalAxes.boxAxes;

  Vec3.add(tmpV, axes.origin, axes.dirC);
  Mat4.targetTo(tmpMat, tmpV, axes.origin, axes.dirB);

  state.currentGroup = 0;

  if (props.variant === 'circle') {
    // the Circle primitive lies in the xz-plane (y as normal), rotate it
    // by 90 degrees around x so it aligns with the plane's local frame
    // (z as normal), matching the orientation used for the rectangle
    Mat4.mul(tmpMat, tmpMat, Mat4.rotX90);
    Mat4.scaleUniformly(tmpMat, tmpMat, props.scaleFactor);
    Mat4.setTranslation(tmpMat, axes.origin);

    const radius = Math.max(Vec3.magnitude(axes.dirA), Vec3.magnitude(axes.dirB));
    const circle = Circle({ radius, segments: props.radialSegments });
    MeshBuilder.addPrimitive(state, tmpMat, circle);
    MeshBuilder.addPrimitiveFlipped(state, tmpMat, circle);
  } else {
    const plane = Plane();
    Mat4.scale(tmpMat, tmpMat, Axes3D.size(tmpV, axes));
    Mat4.scaleUniformly(tmpMat, tmpMat, props.scaleFactor);
    Mat4.setTranslation(tmpMat, axes.origin);

    MeshBuilder.addPrimitive(state, tmpMat, plane);
    MeshBuilder.addPrimitiveFlipped(state, tmpMat, plane);
  }

  return MeshBuilder.getMesh(state);
}

function getPlaneShape(ctx: RuntimeContext, data: PlaneData, props: PlaneProps, shape?: Shape<Mesh>) {
  const mesh = buildPlaneMesh(data, props, shape && shape.geometry);
  const name = getPlaneName(data.locis);
  return ShapeUtils.create(
    name,
    data,
    mesh,
    () => props.color,
    () => 1,
    () => name,
  );
}

//

export type PlaneRepresentation = Representation<PlaneData, PlaneParams>;
export function PlaneRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<PlaneData, PlaneParams>,
): PlaneRepresentation {
  const repr = Representation.createMulti(
    'Plane',
    ctx,
    getParams,
    Representation.StateBuilder,
    PlaneVisuals as unknown as Representation.Def<PlaneData, PlaneParams>,
  );
  repr.setState({ markerActions: MarkerActions.Highlighting });
  return repr;
}
