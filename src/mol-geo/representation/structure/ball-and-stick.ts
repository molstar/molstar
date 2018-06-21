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
import { InterUnitLinkVisual } from './visual/inter-unit-link-cylinder';

export const DefaultBallAndStickProps = {
    ...DefaultElementSphereProps,
    ...DefaultIntraUnitLinkProps,

    sizeTheme: { name: 'uniform', value: 0.25 } as SizeTheme,
}
export type BallAndStickProps = Partial<typeof DefaultBallAndStickProps>

export function BallAndStickRepresentation(): StructureRepresentation<BallAndStickProps> {
    const elmementRepr = StructureRepresentation(ElementSphereVisual)
    const linkRepr = StructureRepresentation(IntraUnitLinkVisual, InterUnitLinkVisual)

    return {
        get renderObjects() {
            // return linkRepr.renderObjects
            return [ ...elmementRepr.renderObjects, ...linkRepr.renderObjects ]
        },
        create: (structure: Structure, props: BallAndStickProps = {} as BallAndStickProps) => {
            const p = Object.assign({}, DefaultBallAndStickProps, props)
            return Task.create('Creating BallAndStickRepresentation', async ctx => {
                await elmementRepr.create(structure, p).runInContext(ctx)
                await linkRepr.create(structure, p).runInContext(ctx)
            })
        },
        update: (props: BallAndStickProps) => {
            return Task.create('Updating BallAndStickRepresentation', async ctx => {
                await elmementRepr.update(props).runInContext(ctx)
                await linkRepr.update(props).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const sphereLoci = elmementRepr.getLoci(pickingId)
            const intraLinkLoci = linkRepr.getLoci(pickingId)
            if (isEmptyLoci(sphereLoci)) {
                return intraLinkLoci
            } else {
                return sphereLoci
            }
        },
        mark: (loci: Loci, action: MarkerAction) => {
            elmementRepr.mark(loci, action)
            linkRepr.mark(loci, action)
        },
        destroy() {
            elmementRepr.destroy()
            linkRepr.destroy()
        }
    }
}