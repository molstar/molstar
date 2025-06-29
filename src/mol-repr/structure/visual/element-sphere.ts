/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsSpheresParams, UnitsVisual, UnitsSpheresVisual, UnitsMeshVisual } from '../units-visual';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { createElementSphereImpostor, ElementIterator, getElementLoci, eachElement, createElementSphereMesh, createStructureElementSphereImpostor, getSerialElementLoci, eachSerialElement, createStructureElementSphereMesh } from './util/element';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Structure } from '../../../mol-model/structure';
import { checkSphereImpostorSupport, StructureGroup } from './util/common';
import { ComplexMeshParams, ComplexMeshVisual, ComplexSpheresParams, ComplexSpheresVisual, ComplexVisual } from '../complex-visual';

export const CommonElementSphereParams = {
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    traceOnly: PD.Boolean(false),
    tryUseImpostor: PD.Boolean(true),
    stride: PD.Numeric(1, { min: 1, max: 100, step: 1 }),
};

//

export const ElementSphereParams = {
    ...UnitsMeshParams,
    ...UnitsSpheresParams,
    ...CommonElementSphereParams,
};
export type ElementSphereParams = typeof ElementSphereParams

export function ElementSphereVisual(materialId: number, structure: Structure, props: PD.Values<ElementSphereParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && checkSphereImpostorSupport(webgl)
        ? ElementSphereImpostorVisual(materialId)
        : ElementSphereMeshVisual(materialId);
}

export function ElementSphereImpostorVisual(materialId: number): UnitsVisual<ElementSphereParams> {
    return UnitsSpheresVisual<ElementSphereParams>({
        defaultProps: PD.getDefaultValues(ElementSphereParams),
        createGeometry: createElementSphereImpostor,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementSphereParams>, currentProps: PD.Values<ElementSphereParams>) => {
            state.createGeometry = (
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.stride !== currentProps.stride
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<ElementSphereParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

export function ElementSphereMeshVisual(materialId: number): UnitsVisual<ElementSphereParams> {
    return UnitsMeshVisual<ElementSphereParams>({
        defaultProps: PD.getDefaultValues(ElementSphereParams),
        createGeometry: createElementSphereMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ElementSphereParams>, currentProps: PD.Values<ElementSphereParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.stride !== currentProps.stride
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<ElementSphereParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}

//

export const StructureElementSphereParams = {
    ...ComplexMeshParams,
    ...ComplexSpheresParams,
    ...CommonElementSphereParams,
};
export type StructureElementSphereParams = typeof ElementSphereParams

export function StructureElementSphereVisual(materialId: number, structure: Structure, props: PD.Values<ElementSphereParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? StructureElementSphereImpostorVisual(materialId)
        : StructureElementSphereMeshVisual(materialId);
}

export function StructureElementSphereImpostorVisual(materialId: number): ComplexVisual<StructureElementSphereParams> {
    return ComplexSpheresVisual<StructureElementSphereParams>({
        defaultProps: PD.getDefaultValues(StructureElementSphereParams),
        createGeometry: createStructureElementSphereImpostor,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureElementSphereParams>, currentProps: PD.Values<StructureElementSphereParams>) => {
            state.createGeometry = (
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.stride !== currentProps.stride
            );
        },
        mustRecreate: (structure: Structure, props: PD.Values<StructureElementSphereParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

export function StructureElementSphereMeshVisual(materialId: number): ComplexVisual<StructureElementSphereParams> {
    return ComplexMeshVisual<StructureElementSphereParams>({
        defaultProps: PD.getDefaultValues(StructureElementSphereParams),
        createGeometry: createStructureElementSphereMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureElementSphereParams>, currentProps: PD.Values<StructureElementSphereParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.stride !== currentProps.stride
            );
        },
        mustRecreate: (structure: Structure, props: PD.Values<StructureElementSphereParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}