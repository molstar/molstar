/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComplexRepresentation, StructureRepresentation, UnitsRepresentation } from '..';
import { ElementSphereVisual, DefaultElementSphereProps } from '../visual/element-sphere';
import { IntraUnitLinkVisual, DefaultIntraUnitLinkProps } from '../visual/intra-unit-link-cylinder';
import { PickingId } from '../../../util/picking';
import { Structure, Unit } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../../util/marker-data';
import { InterUnitLinkVisual } from '../visual/inter-unit-link-cylinder';
import { SizeThemeProps } from 'mol-view/theme/size';

export const DefaultBallAndStickProps = {
    ...DefaultElementSphereProps,
    ...DefaultIntraUnitLinkProps,

    sizeTheme: { name: 'uniform', value: 0.25 } as SizeThemeProps,
    unitKinds: [ Unit.Kind.Atomic ] as Unit.Kind[]
}
export type BallAndStickProps = typeof DefaultBallAndStickProps

export type BallAndStickRepresentation = StructureRepresentation<BallAndStickProps>

export function BallAndStickRepresentation(): BallAndStickRepresentation {
    const elmementRepr = UnitsRepresentation(ElementSphereVisual)
    const intraLinkRepr = UnitsRepresentation(IntraUnitLinkVisual)
    const interLinkRepr = ComplexRepresentation(InterUnitLinkVisual)

    let currentProps: BallAndStickProps
    return {
        get renderObjects() {
            return [ ...elmementRepr.renderObjects, ...intraLinkRepr.renderObjects, ...interLinkRepr.renderObjects ]
        },
        get props() {
            return { ...elmementRepr.props, ...intraLinkRepr.props, ...interLinkRepr.props }
        },
        create: (structure: Structure, props: Partial<BallAndStickProps> = {}) => {
            currentProps = Object.assign({}, DefaultBallAndStickProps, props)
            return Task.create('DistanceRestraintRepresentation', async ctx => {
                await elmementRepr.create(structure, currentProps).runInContext(ctx)
                await intraLinkRepr.create(structure, currentProps).runInContext(ctx)
                await interLinkRepr.create(structure, currentProps).runInContext(ctx)
            })
        },
        update: (props: Partial<BallAndStickProps>) => {
            currentProps = Object.assign(currentProps, props)
            return Task.create('Updating BallAndStickRepresentation', async ctx => {
                await elmementRepr.update(currentProps).runInContext(ctx)
                await intraLinkRepr.update(currentProps).runInContext(ctx)
                await interLinkRepr.update(currentProps).runInContext(ctx)
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