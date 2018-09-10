/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { GaussianSurfaceVisual, DefaultGaussianSurfaceProps } from '../visual/gaussian-surface-mesh';
import { StructureRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { MarkerAction } from '../../../util/marker-data';
import { Loci } from 'mol-model/loci';
import { PickingId } from '../../../util/picking';

export const DefaultSurfaceProps = {
    ...DefaultGaussianSurfaceProps,
}
export type SurfaceProps = typeof DefaultSurfaceProps

export type SurfaceRepresentation = StructureRepresentation<SurfaceProps>

export function SurfaceRepresentation(): SurfaceRepresentation {
    let currentProps: SurfaceProps
    const gaussianRepr = UnitsRepresentation('Gaussian surface', GaussianSurfaceVisual)
    return {
        label: 'Surface',
        get renderObjects() {
            return [ ...gaussianRepr.renderObjects ]
        },
        get props() {
            return { ...gaussianRepr.props }
        },
        createOrUpdate: (props: Partial<SurfaceProps> = {}, structure?: Structure) => {
            currentProps = Object.assign({}, DefaultSurfaceProps, currentProps, props)
            return gaussianRepr.createOrUpdate(currentProps, structure)
        },
        getLoci: (pickingId: PickingId) => {
            return gaussianRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            return gaussianRepr.mark(loci, action)
        },
        destroy() {
            gaussianRepr.destroy()
        }
    }
}