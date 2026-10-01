/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type UnitsVisual, UnitsLinesVisual, UnitsLinesParams } from '../units-visual.js';
import type { VisualContext } from '../../visual.js';
import type { Unit, Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { CommonMolecularSurfaceCalculationParams, computeUnitMolecularSurface, computeStructureMolecularSurface } from './util/molecular-surface.js';
import { computeMarchingCubesLines } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element.js';
import type { VisualUpdateState } from '../../util.js';
import { CommonSurfaceParams } from './util/common.js';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Tensor } from '@molstar/core/math/linear-algebra/tensor';
import { ComplexLinesParams, ComplexLinesVisual, type ComplexVisual } from '../complex-visual.js';

const SharedParams = {
    ...CommonMolecularSurfaceCalculationParams,
    ...CommonSurfaceParams,
    sizeFactor: PD.Numeric(1.5, { min: 0, max: 10, step: 0.1 }),
};
type SharedParams = typeof SharedParams

export const MolecularSurfaceWireframeParams = {
    ...UnitsLinesParams,
    ...SharedParams,
};
export type MolecularSurfaceWireframeParams = typeof MolecularSurfaceWireframeParams
export type MolecularSurfaceWireframeProps = PD.Values<MolecularSurfaceWireframeParams>

export const StructureMolecularSurfaceWireframeParams = {
    ...ComplexLinesParams,
    ...SharedParams,
};
export type StructureMolecularSurfaceWireframeParams = typeof StructureMolecularSurfaceWireframeParams
export type StructureMolecularSurfaceWireframeProps = PD.Values<StructureMolecularSurfaceWireframeParams>

//

async function createMolecularSurfaceWireframe(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: MolecularSurfaceWireframeProps, lines?: Lines): Promise<Lines> {
    const { transform, field, idField, maxRadius } = await computeUnitMolecularSurface(structure, unit, theme.size, props).runInContext(ctx.runtime);

    const params = {
        isoLevel: props.probeRadius,
        scalarField: props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export function MolecularSurfaceWireframeVisual(materialId: number): UnitsVisual<MolecularSurfaceWireframeParams> {
    return UnitsLinesVisual<MolecularSurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(MolecularSurfaceWireframeParams),
        createGeometry: createMolecularSurfaceWireframe,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<MolecularSurfaceWireframeParams>, currentProps: PD.Values<MolecularSurfaceWireframeParams>) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.probeRadius !== currentProps.probeRadius ||
                newProps.probePositions !== currentProps.probePositions ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill
            );
        }
    }, materialId);
}

//

async function createStructureMolecularSurfaceWireframe(ctx: VisualContext, structure: Structure, theme: Theme, props: StructureMolecularSurfaceWireframeProps, lines?: Lines): Promise<Lines> {
    const { transform, field, idField, maxRadius } = await computeStructureMolecularSurface(structure, theme.size, props).runInContext(ctx.runtime);

    const params = {
        isoLevel: props.probeRadius,
        scalarField: props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
        idField
    };
    const wireframe = await computeMarchingCubesLines(params, lines).runAsChild(ctx.runtime);

    Lines.transform(wireframe, transform);

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius);
    wireframe.setBoundingSphere(sphere);

    return wireframe;
}

export function StructureMolecularSurfaceWireframeVisual(materialId: number): ComplexVisual<StructureMolecularSurfaceWireframeParams> {
    return ComplexLinesVisual<StructureMolecularSurfaceWireframeParams>({
        defaultProps: PD.getDefaultValues(StructureMolecularSurfaceWireframeParams),
        createGeometry: createStructureMolecularSurfaceWireframe,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureMolecularSurfaceWireframeParams>, currentProps: PD.Values<StructureMolecularSurfaceWireframeParams>) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.probeRadius !== currentProps.probeRadius ||
                newProps.probePositions !== currentProps.probePositions ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill
            );
        }
    }, materialId);
}
