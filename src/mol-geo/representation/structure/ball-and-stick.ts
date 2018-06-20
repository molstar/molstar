/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation } from '.';
import { ElementSphereVisual, DefaultElementSphereProps } from './visual/element-sphere';
import { IntraUnitLinkVisual, DefaultIntraUnitLinkProps } from './visual/intra-unit-link-cylinder';
import { PickingId } from '../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';
import { SizeTheme } from '../../theme';

export const DefaultBallAndStickProps = {
    ...DefaultElementSphereProps,
    ...DefaultIntraUnitLinkProps,

    sizeTheme: { name: 'uniform', value: 0.25 } as SizeTheme,
}
export type BallAndStickProps = Partial<typeof DefaultBallAndStickProps>

export function BallAndStickRepresentation(): StructureRepresentation<BallAndStickProps> {
    const sphereRepr = StructureRepresentation(ElementSphereVisual)
    const intraLinkRepr = StructureRepresentation(IntraUnitLinkVisual)

    return {
        get renderObjects() {
            return [ ...sphereRepr.renderObjects, ...intraLinkRepr.renderObjects ]
        },
        create: (structure: Structure, props: BallAndStickProps = {} as BallAndStickProps) => {
            const p = Object.assign({}, DefaultBallAndStickProps, props)
            return Task.create('Creating BallAndStickRepresentation', async ctx => {
                await sphereRepr.create(structure, p).runInContext(ctx)
                await intraLinkRepr.create(structure, p).runInContext(ctx)
            })
        },
        update: (props: BallAndStickProps) => {
            return Task.create('Updating BallAndStickRepresentation', async ctx => {
                await sphereRepr.update(props).runInContext(ctx)
                await intraLinkRepr.update(props).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const sphereLoci = sphereRepr.getLoci(pickingId)
            const intraLinkLoci = intraLinkRepr.getLoci(pickingId)
            if (isEmptyLoci(sphereLoci)) {
                return intraLinkLoci
            } else {
                return sphereLoci
            }
        },
        mark: (loci: Loci, action: MarkerAction) => {
            sphereRepr.mark(loci, action)
            intraLinkRepr.mark(loci, action)
        },
        destroy() {
            sphereRepr.destroy()
            intraLinkRepr.destroy()
        }
    }
}