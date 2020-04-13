/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { CylinderProps } from '../../../mol-geo/primitive/cylinder';
import { PolymerBackboneIterator } from './util/polymer';
import { OrderedSet } from '../../../mol-data/int';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { ElementIterator, getElementLoci, eachElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';

export const PolymerBackboneCylinderParams = {
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
};
export const DefaultPolymerBackboneCylinderProps = PD.getDefaultValues(PolymerBackboneCylinderParams);
export type PolymerBackboneCylinderProps = typeof DefaultPolymerBackboneCylinderProps

// TODO do group id based on polymer index not element index
function createPolymerBackboneCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerBackboneCylinderProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length;
    if (!polymerElementCount) return Mesh.createEmpty(mesh);

    const { radialSegments, sizeFactor } = props;

    const vertexCountEstimate = radialSegments * 2 * polymerElementCount * 2;
    const builderState = MeshBuilder.createState(vertexCountEstimate, vertexCountEstimate / 10, mesh);

    const { elements } = unit;
    const pos = unit.conformation.invariantPosition;
    const pA = Vec3.zero();
    const pB = Vec3.zero();
    const cylinderProps: CylinderProps = { radiusTop: 1, radiusBottom: 1, radialSegments };

    const polymerBackboneIt = PolymerBackboneIterator(structure, unit);
    while (polymerBackboneIt.hasNext) {
        const { centerA, centerB } = polymerBackboneIt.move();
        pos(centerA.element, pA);
        pos(centerB.element, pB);

        cylinderProps.radiusTop = cylinderProps.radiusBottom = theme.size.size(centerA) * sizeFactor;
        builderState.currentGroup = OrderedSet.indexOf(elements, centerA.element);
        addCylinder(builderState, pA, pB, 0.5, cylinderProps);

        cylinderProps.radiusTop = cylinderProps.radiusBottom = theme.size.size(centerB) * sizeFactor;
        builderState.currentGroup = OrderedSet.indexOf(elements, centerB.element);
        addCylinder(builderState, pB, pA, 0.5, cylinderProps);
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export const PolymerBackboneParams = {
    ...UnitsMeshParams,
    ...PolymerBackboneCylinderParams,
};
export type PolymerBackboneParams = typeof PolymerBackboneParams

export function PolymerBackboneVisual(materialId: number): UnitsVisual<PolymerBackboneParams> {
    return UnitsMeshVisual<PolymerBackboneParams>({
        defaultProps: PD.getDefaultValues(PolymerBackboneParams),
        createGeometry: createPolymerBackboneCylinderMesh,
        // TODO create a specialized location iterator
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerBackboneParams>, currentProps: PD.Values<PolymerBackboneParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        }
    }, materialId);
}