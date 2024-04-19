/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure, ElementIndex, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { eachPolymerElement, getPolymerElementLoci, PolymerLocationIterator } from './util/polymer';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, UnitsSpheresVisual, UnitsSpheresParams } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { sphereVertexCount } from '../../../mol-geo/primitive/sphere';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../../mol-geo/geometry/spheres/spheres-builder';
import { eachPolymerBackboneElement } from './util/polymer/backbone';
import { StructureGroup } from './util/common';

export const PolymerBackboneSphereParams = {
    ...UnitsMeshParams,
    ...UnitsSpheresParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true),
};
export type PolymerBackboneSphereParams = typeof PolymerBackboneSphereParams

export function PolymerBackboneSphereVisual(materialId: number, structure: Structure, props: PD.Values<PolymerBackboneSphereParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? PolymerBackboneSphereImpostorVisual(materialId)
        : PolymerBackboneSphereMeshVisual(materialId);
}

interface PolymerBackboneSphereProps {
    detail: number,
    sizeFactor: number,
}

function createPolymerBackboneSphereImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerBackboneSphereProps, spheres?: Spheres) {
    const polymerElementCount = unit.polymerElements.length;
    if (!polymerElementCount) return Spheres.createEmpty(spheres);

    const builder = SpheresBuilder.create(polymerElementCount, polymerElementCount / 2, spheres);

    const c = unit.conformation;
    const p = Vec3();

    const add = (index: ElementIndex, group: number) => {
        c.invariantPosition(index, p);
        builder.add(p[0], p[1], p[2], group);
    };

    eachPolymerBackboneElement(unit, add);

    const s = builder.getSpheres();

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    s.setBoundingSphere(sphere);

    return s;
}

export function PolymerBackboneSphereImpostorVisual(materialId: number): UnitsVisual<PolymerBackboneSphereParams> {
    return UnitsSpheresVisual<PolymerBackboneSphereParams>({
        defaultProps: PD.getDefaultValues(PolymerBackboneSphereParams),
        createGeometry: createPolymerBackboneSphereImpostor,
        createLocationIterator: (structureGroup: StructureGroup) => PolymerLocationIterator.fromGroup(structureGroup),
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerBackboneSphereParams>, currentProps: PD.Values<PolymerBackboneSphereParams>) => { },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<PolymerBackboneSphereParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

function createPolymerBackboneSphereMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PolymerBackboneSphereProps, mesh?: Mesh) {
    const polymerElementCount = unit.polymerElements.length;
    if (!polymerElementCount) return Mesh.createEmpty(mesh);

    const { detail, sizeFactor } = props;

    const vertexCount = polymerElementCount * sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);

    const c = unit.conformation;
    const p = Vec3();
    const center = StructureElement.Location.create(structure, unit);

    const add = (index: ElementIndex, group: number) => {
        center.element = index;
        c.invariantPosition(center.element, p);
        builderState.currentGroup = group;
        addSphere(builderState, p, theme.size.size(center) * sizeFactor, detail);
    };

    eachPolymerBackboneElement(unit, add);

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export function PolymerBackboneSphereMeshVisual(materialId: number): UnitsVisual<PolymerBackboneSphereParams> {
    return UnitsMeshVisual<PolymerBackboneSphereParams>({
        defaultProps: PD.getDefaultValues(PolymerBackboneSphereParams),
        createGeometry: createPolymerBackboneSphereMesh,
        createLocationIterator: (structureGroup: StructureGroup) => PolymerLocationIterator.fromGroup(structureGroup),
        getLoci: getPolymerElementLoci,
        eachLocation: eachPolymerElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PolymerBackboneSphereParams>, currentProps: PD.Values<PolymerBackboneSphereParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<PolymerBackboneSphereParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}