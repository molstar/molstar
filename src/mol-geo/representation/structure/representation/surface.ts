/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { GaussianSurfaceVisual, DefaultGaussianSurfaceProps } from '../visual/gaussian-surface-mesh';
import { StructureRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { MarkerAction } from '../../../geometry/marker-data';
import { Loci } from 'mol-model/loci';
import { PickingId } from '../../../geometry/picking';
import { Task } from 'mol-task';
import { DefaultGaussianWireframeProps, GaussianWireframeVisual } from '../visual/gaussian-surface-wireframe';
import { ColorThemeProps } from 'mol-view/theme/color';

export const DefaultSurfaceProps = {
    ...DefaultGaussianSurfaceProps,
    ...DefaultGaussianWireframeProps,

    visuals: { surface: true, wireframe: false },
    colorTheme: { name: 'element-symbol' } as ColorThemeProps,
}
export type SurfaceProps = typeof DefaultSurfaceProps

export type SurfaceRepresentation = StructureRepresentation<SurfaceProps>

export function SurfaceRepresentation(): SurfaceRepresentation {
    let currentProps: SurfaceProps
    let currentStructure: Structure
    const gaussianSurfaceRepr = UnitsRepresentation('Gaussian surface', GaussianSurfaceVisual)
    const gaussianWireframeRepr = UnitsRepresentation('Gaussian wireframe', GaussianWireframeVisual)
    return {
        label: 'Surface',
        get renderObjects() {
            const renderObjects = []
            if (currentProps.visuals.surface) renderObjects.push(...gaussianSurfaceRepr.renderObjects)
            if (currentProps.visuals.wireframe) renderObjects.push(...gaussianWireframeRepr.renderObjects)
            return renderObjects
        },
        get props() {
            return { ...gaussianSurfaceRepr.props, ...gaussianWireframeRepr.props, visuals: currentProps.visuals }
        },
        createOrUpdate: (props: Partial<SurfaceProps> = {}, structure?: Structure) => {
            currentProps = Object.assign({}, DefaultSurfaceProps, currentProps, props)
            if (structure) currentStructure = structure
            return Task.create('Creating SurfaceRepresentation', async ctx => {
                if (currentProps.visuals.surface) await gaussianSurfaceRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
                if (currentProps.visuals.wireframe) await gaussianWireframeRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            return gaussianSurfaceRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            return gaussianSurfaceRepr.mark(loci, action)
        },
        destroy() {
            gaussianSurfaceRepr.destroy()
            gaussianWireframeRepr.destroy()
        }
    }
}