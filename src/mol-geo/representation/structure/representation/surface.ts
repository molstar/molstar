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
import { GaussianDensityPointVisual, DefaultGaussianDensityPointProps } from '../visual/gaussian-density-point';

export const DefaultSurfaceProps = {
    ...DefaultGaussianSurfaceProps,
    ...DefaultGaussianDensityPointProps,
}
export type SurfaceProps = typeof DefaultSurfaceProps

export type SurfaceRepresentation = StructureRepresentation<SurfaceProps>

export function SurfaceRepresentation(): SurfaceRepresentation {
    let currentProps: SurfaceProps
    const gaussianSurfaceRepr = UnitsRepresentation('Gaussian surface', GaussianSurfaceVisual)
    const gaussianPointRepr = UnitsRepresentation('Gaussian point grid', GaussianDensityPointVisual)
    return {
        label: 'Surface',
        get renderObjects() {
            return [ ...gaussianSurfaceRepr.renderObjects, ...gaussianPointRepr.renderObjects ]
        },
        get props() {
            return { ...gaussianSurfaceRepr.props, ...gaussianPointRepr.props }
        },
        createOrUpdate: (props: Partial<SurfaceProps> = {}, structure?: Structure) => {
            currentProps = Object.assign({}, DefaultSurfaceProps, currentProps, props)
            return Task.create('Creating SurfaceRepresentation', async ctx => {
                await gaussianSurfaceRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
                // await gaussianPointRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
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
            gaussianPointRepr.destroy()
        }
    }
}