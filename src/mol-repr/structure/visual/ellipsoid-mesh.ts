/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../../../mol-repr/structure/units-visual';
import { ElementIterator, getElementLoci, eachElement } from '../../../mol-repr/structure/visual/util/element';
import { VisualUpdateState } from '../../../mol-repr/util';
import { VisualContext } from '../../../mol-repr/visual';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { sphereVertexCount } from '../../../mol-geo/primitive/sphere';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3, Mat3, Tensor, EPSILON } from '../../../mol-math/linear-algebra';
import { isHydrogen } from '../../../mol-repr/structure/visual/util/common';
import { addEllipsoid } from '../../../mol-geo/geometry/mesh/builder/ellipsoid';
import { AtomSiteAnisotrop } from '../../../mol-model-formats/structure/property/anisotropic';
import { equalEps } from '../../../mol-math/linear-algebra/3d/common';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { Sphere3D } from '../../../mol-math/geometry';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

export const EllipsoidMeshParams = {
    ...UnitsMeshParams,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    ignoreHydrogens: PD.Boolean(false),
};
export type EllipsoidMeshParams = typeof EllipsoidMeshParams

export function EllipsoidMeshVisual(materialId: number): UnitsVisual<EllipsoidMeshParams> {
    return UnitsMeshVisual<EllipsoidMeshParams>({
        defaultProps: PD.getDefaultValues(EllipsoidMeshParams),
        createGeometry: createEllipsoidMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<EllipsoidMeshParams>, currentProps: PD.Values<EllipsoidMeshParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens
            );
        }
    }, materialId);
}

//

export interface EllipsoidMeshProps {
    detail: number,
    sizeFactor: number,
    ignoreHydrogens: boolean,
}

export function createEllipsoidMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: EllipsoidMeshProps, mesh?: Mesh): Mesh {
    const { detail, sizeFactor } = props;

    const { elements, model } = unit;
    const elementCount = elements.length;
    const vertexCount = elementCount * sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);

    const atomSiteAnisotrop = AtomSiteAnisotrop.Provider.get(model);
    if (!atomSiteAnisotrop) return Mesh.createEmpty(mesh);

    const v = Vec3();
    const mat = Mat3();
    const eigvals = Vec3();
    const eigvec1 = Vec3();
    const eigvec2 = Vec3();
    const { elementToAnsiotrop, data } = atomSiteAnisotrop;
    const { U } = data;
    const space = data._schema.U.space;
    const pos = unit.conformation.invariantPosition;
    const l = StructureElement.Location.create(structure);
    l.unit = unit;

    for (let i = 0; i < elementCount; i++) {
        const ei = elements[i];
        const ai = elementToAnsiotrop[ei];
        if (ai === -1) continue;
        if (props.ignoreHydrogens && isHydrogen(unit, ei)) continue;

        l.element = ei;
        pos(ei, v);

        builderState.currentGroup = i;
        Tensor.toMat3(mat, space, U.value(ai));
        Mat3.symmtricFromLower(mat, mat);
        Mat3.symmetricEigenvalues(eigvals, mat);
        Mat3.eigenvector(eigvec1, mat, eigvals[1]);
        Mat3.eigenvector(eigvec2, mat, eigvals[2]);
        for (let j = 0; j < 3; ++j) {
            // show 50% probability surface, needs sqrt as U matrix is in angstrom-squared
            // take abs of eigenvalue to avoid reflection
            // multiply by given size-factor
            eigvals[j] = sizeFactor * 1.5958 * Math.sqrt(Math.abs(eigvals[j]));
        }

        if (equalEps(eigvals[0], eigvals[1], EPSILON) && equalEps(eigvals[1], eigvals[2], EPSILON)) {
            addSphere(builderState, v, eigvals[0], detail);
        } else {
            addEllipsoid(builderState, v, eigvec2, eigvec1, eigvals, detail);
        }
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}