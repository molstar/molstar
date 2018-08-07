/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation, StructureUnitsRepresentation } from '.';
import { ElementSphereVisual, DefaultElementSphereProps } from './visual/element-sphere';
import { IntraUnitLinkVisual, DefaultIntraUnitLinkProps } from './visual/intra-unit-link-cylinder';
import { PickingId } from '../../util/picking';
import { Structure, Unit } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';
import { SizeTheme } from '../../theme';
import { InterUnitLinkVisual } from './visual/inter-unit-link-cylinder';

export const DefaultBallAndStickProps = {
    ...DefaultElementSphereProps,
    ...DefaultIntraUnitLinkProps,

    sizeTheme: { name: 'uniform', value: 0.25 } as SizeTheme,
    unitKinds: [ Unit.Kind.Atomic ] as Unit.Kind[]
}
export type BallAndStickProps = Partial<typeof DefaultBallAndStickProps>

export function BallAndStickRepresentation(): StructureRepresentation<BallAndStickProps> {
    const elmementRepr = StructureUnitsRepresentation(ElementSphereVisual)
    const intraLinkRepr = StructureUnitsRepresentation(IntraUnitLinkVisual)
    const interLinkRepr = StructureRepresentation(InterUnitLinkVisual)

    return {
        get renderObjects() {
            return [ ...elmementRepr.renderObjects, ...intraLinkRepr.renderObjects, ...interLinkRepr.renderObjects ]
        },
        get props() {
            return { ...elmementRepr.props, ...intraLinkRepr.props, ...interLinkRepr.props }
        },
        create: (structure: Structure, props: BallAndStickProps = {} as BallAndStickProps) => {
            const p = Object.assign({}, DefaultBallAndStickProps, props)
            return Task.create('DistanceRestraintRepresentation', async ctx => {
                await elmementRepr.create(structure, p).runInContext(ctx)
                await intraLinkRepr.create(structure, p).runInContext(ctx)
                await interLinkRepr.create(structure, p).runInContext(ctx)
            })
        },
        update: (props: BallAndStickProps) => {
            const p = Object.assign({}, props)
            return Task.create('Updating BallAndStickRepresentation', async ctx => {
                await elmementRepr.update(p).runInContext(ctx)
                await intraLinkRepr.update(p).runInContext(ctx)
                await interLinkRepr.update(p).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const sphereLoci = elmementRepr.getLoci(pickingId)
            const intraLinkLoci = intraLinkRepr.getLoci(pickingId)
            const interLinkLoci = interLinkRepr.getLoci(pickingId)
            if (isEmptyLoci(sphereLoci)) {
                if (isEmptyLoci(intraLinkLoci)) {
                    return interLinkLoci
                } else {
                    return intraLinkLoci
                }
            } else {
                return sphereLoci
            }
        },
        mark: (loci: Loci, action: MarkerAction) => {
            elmementRepr.mark(loci, action)
            intraLinkRepr.mark(loci, action)
            interLinkRepr.mark(loci, action)
        },
        destroy() {
            elmementRepr.destroy()
            intraLinkRepr.destroy()
            interLinkRepr.destroy()
        }
    }
}