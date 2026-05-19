/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Prihoda <david.prihoda@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { OctagonalPyramid } from '../../../mol-geo/primitive/pyramid';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { createCurveSegmentState, PolymerTraceIterator, interpolateCurveSegment, PolymerLocationIterator, getPolymerElementLoci, eachPolymerElement } from './util/polymer';
import { isNucleic } from '../../../mol-model/structure/model/types';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';
import { StructureGroup } from './util/common';

// In PolymerTraceIterator the control points are shifted by one residue:
//   v.p2 = current residue position
//   v.p1 = previous residue position
//   v.p3 = next residue position

const t = Mat4.identity();
const sVec = Vec3.zero();
const offsetVec = Vec3.zero();
const n0 = Vec3.zero();
const n1 = Vec3.zero();
const upVec = Vec3.zero();

const cone = OctagonalPyramid();

type TerminusType = 'nterm' | 'cterm'

function createPolymerTerminusConeMesh(terminusType: TerminusType, ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<PolymerTerminaParams>, mesh?: Mesh) {
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

        const shouldDraw = terminusType === 'nterm' ? v.first : v.last;
        if (!shouldDraw) {
            ++i;
            continue;
        }

        builderState.currentGroup = i;

        const isNucleicType = isNucleic(v.moleculeType);
        const size = theme.size.size(v.center) * sizeFactor;
        const height = 2.0 * size;
        const radius = 1.2 * size;

        if (terminusType === 'nterm') {
            // v.p2 = current (first) residue, v.p3 = next residue.
            // At segment start the interpolated normals are degenerate — compute upVec fresh.
            Vec3.sub(n0, v.p3, v.p2);
            const chainLen = Vec3.magnitude(n0);
            Vec3.set(upVec, 0, 1, 0);
            if (chainLen > 0.01) {
                Vec3.scale(n0, n0, 1 / chainLen);
                if (Math.abs(Vec3.dot(n0, upVec)) > 0.9) Vec3.set(upVec, 1, 0, 0);
            }
            // Z = normalize(v.p3 - v.p2) = forward; tip (+Z) points into chain.
            // Base (-Z) extends outward past chain start.
            Mat4.targetTo(t, v.p3, v.p2, upVec);
        } else {
            // v.p2 = current (last) residue, v.p1 = previous residue.
            const tension = isNucleicType ? 0.5 : 0.9;
            const shift = isNucleicType ? 0.3 : 0.5;
            interpolateCurveSegment(state, v, tension, shift);
            const vectors = isNucleicType ? binormalVectors : normalVectors;
            Vec3.fromArray(n0, vectors, 0);
            Vec3.fromArray(n1, vectors, 3);
            Vec3.normalize(upVec, Vec3.add(upVec, n0, n1));
            // Z = normalize(v.p2 - v.p1) = forward; tip (+Z) points away from chain end.
            Mat4.targetTo(t, v.p2, v.p1, upVec);
        }

        // For nterm: Z is forward, base at -Z. Shift center by -height/2 along Z so base sits at v.p2.
        // For cterm: Z is forward, tip at +Z. Shift center by +height/2 along Z so base sits at v.p2.
        const zSign = terminusType === 'nterm' ? -0.5 : 0.5;
        Vec3.set(offsetVec,
            v.p2[0] + t[8] * height * zSign,
            v.p2[1] + t[9] * height * zSign,
            v.p2[2] + t[10] * height * zSign
        );
        Mat4.scale(t, t, Vec3.set(sVec, radius, radius, height));
        Mat4.setTranslation(t, offsetVec);
        MeshBuilder.addPrimitive(builderState, t, cone);

        ++i;
    }

    const m = MeshBuilder.getMesh(builderState);
    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export const PolymerTerminaParams = {
    ...UnitsMeshParams,
    sizeFactor: PD.Numeric(1.5, { min: 0.01, max: 10, step: 0.01 }),
};
export type PolymerTerminaParams = typeof PolymerTerminaParams

export function PolymerNtermConeVisual(materialId: number): UnitsVisual<PolymerTerminaParams> {
    return UnitsMeshVisual<PolymerTerminaParams>({
        defaultProps: PD.getDefaultValues(PolymerTerminaParams),
        createGeometry: createPolymerTerminusConeMesh.bind(null, 'nterm'),
        createLocationIterator: (structureGroup: StructureGroup) => PolymerLocationIterator.fromGroup(structureGroup),
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerTerminaParams>, currentProps: PD.Values<PolymerTerminaParams>) => {
            state.createGeometry = newProps.sizeFactor !== currentProps.sizeFactor;
        }
    }, materialId);
}

export function PolymerCtermConeVisual(materialId: number): UnitsVisual<PolymerTerminaParams> {
    return UnitsMeshVisual<PolymerTerminaParams>({
        defaultProps: PD.getDefaultValues(PolymerTerminaParams),
        createGeometry: createPolymerTerminusConeMesh.bind(null, 'cterm'),
        createLocationIterator: (structureGroup: StructureGroup) => PolymerLocationIterator.fromGroup(structureGroup),
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerTerminaParams>, currentProps: PD.Values<PolymerTerminaParams>) => {
            state.createGeometry = newProps.sizeFactor !== currentProps.sizeFactor;
        }
    }, materialId);
}
