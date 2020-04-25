/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsSpheresParams, UnitsVisual, UnitsSpheresVisual, UnitsMeshVisual } from '../units-visual';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { createElementSphereImpostor, ElementIterator, getElementLoci, eachElement, createElementSphereMesh } from './util/element';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

export const ElementSphereParams = {
    ...UnitsMeshParams,
    ...UnitsSpheresParams,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    ignoreHydrogens: PD.Boolean(false),
    traceOnly: PD.Boolean(false),
};
export type ElementSphereParams = typeof ElementSphereParams

export function getElementSphereVisual(webgl?: WebGLContext) {
    return webgl && webgl.extensions.fragDepth ? ElementSphereImpostorVisual : ElementSphereMeshVisual;
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
                newProps.traceOnly !== currentProps.traceOnly
            );
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
                newProps.traceOnly !== currentProps.traceOnly
            );
        }
    }, materialId);
}