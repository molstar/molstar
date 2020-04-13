/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Wedge } from '../../../mol-geo/primitive/wedge';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { createCurveSegmentState, PolymerTraceIterator, interpolateCurveSegment, PolymerLocationIterator, getPolymerElementLoci, eachPolymerElement } from './util/polymer';
import { isNucleic, SecondaryStructureType } from '../../../mol-model/structure/model/types';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';

const t = Mat4.identity();
const sVec = Vec3.zero();
const n0 = Vec3.zero();
const n1 = Vec3.zero();
const upVec = Vec3.zero();

const depthFactor = 4;
const widthFactor = 4;
const heightFactor = 6;

const wedge = Wedge();

export const PolymerDirectionWedgeParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
};
export const DefaultPolymerDirectionWedgeProps = PD.getDefaultValues(PolymerDirectionWedgeParams);
export type PolymerDirectionWedgeProps = typeof DefaultPolymerDirectionWedgeProps

function createPolymerDirectionWedgeMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerDirectionWedgeProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length;

    if (!polymerElementCount) return Mesh.createEmpty(mesh);
    const { sizeFactor } = props;

    const vertexCount = polymerElementCount * 24;
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 10, mesh);
    const linearSegments = 1;

    const state = createCurveSegmentState(linearSegments);
    const { normalVectors, binormalVectors } = state;

    let i = 0;
    const polymerTraceIt = PolymerTraceIterator(unit, structure);
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move();
        builderState.currentGroup = i;

        const isNucleicType = isNucleic(v.moleculeType);
        const isSheet = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta);
        const tension = (isNucleicType || isSheet) ? 0.5 : 0.9;
        const shift = isNucleicType ? 0.3 : 0.5;

        interpolateCurveSegment(state, v, tension, shift);

        if ((isSheet && !v.secStrucLast) || !isSheet) {
            const size = theme.size.size(v.center) * sizeFactor;
            const depth = depthFactor * size;
            const width = widthFactor * size;
            const height = heightFactor * size;

            const vectors = isNucleicType ? binormalVectors : normalVectors;
            Vec3.fromArray(n0, vectors, 0);
            Vec3.fromArray(n1, vectors, 3);
            Vec3.normalize(upVec, Vec3.add(upVec, n0, n1));

            Mat4.targetTo(t, v.p3, v.p1, upVec);
            Mat4.mul(t, t, Mat4.rotY90Z180);
            Mat4.scale(t, t, Vec3.set(sVec, height, width, depth));
            Mat4.setTranslation(t, v.p2);
            MeshBuilder.addPrimitive(builderState, t, wedge);
        }

        ++i;
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export const PolymerDirectionParams = {
    ...UnitsMeshParams,
    ...PolymerDirectionWedgeParams
};
export type PolymerDirectionParams = typeof PolymerDirectionParams

export function PolymerDirectionVisual(materialId: number): UnitsVisual<PolymerDirectionParams> {
    return UnitsMeshVisual<PolymerDirectionParams>({
        defaultProps: PD.getDefaultValues(PolymerDirectionParams),
        createGeometry: createPolymerDirectionWedgeMesh,
        createLocationIterator: PolymerLocationIterator.fromGroup,
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerDirectionParams>, currentProps: PD.Values<PolymerDirectionParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor
            );
        }
    }, materialId);
}