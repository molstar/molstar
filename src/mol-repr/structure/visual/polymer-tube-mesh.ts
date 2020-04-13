/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { createCurveSegmentState, PolymerTraceIterator, interpolateCurveSegment, interpolateSizes, PolymerLocationIterator, getPolymerElementLoci, eachPolymerElement, HelixTension, StandardTension, StandardShift, NucleicShift, OverhangFactor } from './util/polymer';
import { isNucleic, SecondaryStructureType } from '../../../mol-model/structure/model/types';
import { addTube } from '../../../mol-geo/geometry/mesh/builder/tube';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { addSheet } from '../../../mol-geo/geometry/mesh/builder/sheet';
import { addRibbon } from '../../../mol-geo/geometry/mesh/builder/ribbon';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';

export const PolymerTubeMeshParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    linearSegments: PD.Numeric(8, { min: 1, max: 48, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
};
export const DefaultPolymerTubeMeshProps = PD.getDefaultValues(PolymerTubeMeshParams);
export type PolymerTubeMeshProps = typeof DefaultPolymerTubeMeshProps

const tmpV1 = Vec3();

function createPolymerTubeMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerTubeMeshProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length;

    if (!polymerElementCount) return Mesh.createEmpty(mesh);
    const { sizeFactor, detail, linearSegments, radialSegments } = props;

    const vertexCount = linearSegments * radialSegments * polymerElementCount + (radialSegments + 1) * polymerElementCount * 2;
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 10, mesh);

    const state = createCurveSegmentState(linearSegments);
    const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues } = state;

    let i = 0;
    const polymerTraceIt = PolymerTraceIterator(unit, structure);
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move();
        builderState.currentGroup = i;

        const isNucleicType = isNucleic(v.moleculeType);
        const isHelix = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Helix);
        const tension = isHelix ? HelixTension : StandardTension;
        const shift = isNucleicType ? NucleicShift : StandardShift;

        interpolateCurveSegment(state, v, tension, shift);

        const startCap = v.secStrucFirst || v.coarseBackboneFirst || v.first;
        const endCap = v.secStrucLast || v.coarseBackboneLast || v.last;

        let s0 = theme.size.size(v.centerPrev) * sizeFactor;
        let s1 = theme.size.size(v.center) * sizeFactor;
        let s2 = theme.size.size(v.centerNext) * sizeFactor;

        interpolateSizes(state, s0, s1, s2, s0, s1, s2, shift);

        let segmentCount = linearSegments;
        if (v.initial) {
            segmentCount = Math.max(Math.round(linearSegments * shift), 1);
            const offset = linearSegments - segmentCount;
            curvePoints.copyWithin(0, offset * 3);
            binormalVectors.copyWithin(0, offset * 3);
            normalVectors.copyWithin(0, offset * 3);
            widthValues.copyWithin(0, offset * 3);
            heightValues.copyWithin(0, offset * 3);
            Vec3.fromArray(tmpV1, curvePoints, 3);
            Vec3.normalize(tmpV1, Vec3.sub(tmpV1, v.p2, tmpV1));
            Vec3.scaleAndAdd(tmpV1, v.p2, tmpV1, s1 * OverhangFactor);
            Vec3.toArray(tmpV1, curvePoints, 0);
        } else if (v.final) {
            segmentCount = Math.max(Math.round(linearSegments * (1 - shift)), 1);
            Vec3.fromArray(tmpV1, curvePoints, segmentCount * 3 - 3);
            Vec3.normalize(tmpV1, Vec3.sub(tmpV1, v.p2, tmpV1));
            Vec3.scaleAndAdd(tmpV1, v.p2, tmpV1, s1 * OverhangFactor);
            Vec3.toArray(tmpV1, curvePoints, segmentCount * 3);
        }

        if (v.initial === true && v.final === true) {
            addSphere(builderState, v.p2, s1 * 2, detail);
        } else if (radialSegments === 2) {
            addRibbon(builderState, curvePoints, normalVectors, binormalVectors, segmentCount, widthValues, heightValues, 0);
        } else if (radialSegments === 4) {
            addSheet(builderState, curvePoints, normalVectors, binormalVectors, segmentCount, widthValues, heightValues, 0, startCap, endCap);
        } else {
            addTube(builderState, curvePoints, normalVectors, binormalVectors, segmentCount, radialSegments, widthValues, heightValues, 1, startCap, endCap);
        }

        ++i;
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export const PolymerTubeParams = {
    ...UnitsMeshParams,
    ...PolymerTubeMeshParams
};
export type PolymerTubeParams = typeof PolymerTubeParams

export function PolymerTubeVisual(materialId: number): UnitsVisual<PolymerTubeParams> {
    return UnitsMeshVisual<PolymerTubeParams>({
        defaultProps: PD.getDefaultValues(PolymerTubeParams),
        createGeometry: createPolymerTubeMesh,
        createLocationIterator: PolymerLocationIterator.fromGroup,
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerTubeParams>, currentProps: PD.Values<PolymerTubeParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail ||
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        }
    }, materialId);
}