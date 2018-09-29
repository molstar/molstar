/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { GaussianSurfaceVisual, GaussianSurfaceParams } from '../visual/gaussian-surface-mesh';
import { StructureRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { MarkerAction } from '../../../geometry/marker-data';
import { Loci } from 'mol-model/loci';
import { PickingId } from '../../../geometry/picking';
import { Task } from 'mol-task';
import { GaussianWireframeVisual, GaussianWireframeParams } from '../visual/gaussian-surface-wireframe';
import { getQualityProps } from '../../util';
import { paramDefaultValues, MultiSelectParam } from 'mol-view/parameter';

const VisualOptions = [['surface', 'Surface'], ['wireframe', 'Wireframe']] as [string, string][]

export const MolecularSurfaceParams = {
    ...GaussianSurfaceParams,
    ...GaussianWireframeParams,

    visuals: MultiSelectParam<string>('Visuals', '', ['surface'], VisualOptions)
}
export const DefaultMolecularSurfaceProps = paramDefaultValues(MolecularSurfaceParams)
export type MolecularSurfaceProps = typeof DefaultMolecularSurfaceProps

export type MolecularSurfaceRepresentation = StructureRepresentation<MolecularSurfaceProps>

export function MolecularSurfaceRepresentation(): MolecularSurfaceRepresentation {
    let currentProps: MolecularSurfaceProps
    let currentStructure: Structure
    const gaussianSurfaceRepr = UnitsRepresentation('Gaussian surface', GaussianSurfaceVisual)
    const gaussianWireframeRepr = UnitsRepresentation('Gaussian wireframe', GaussianWireframeVisual)
    return {
        label: 'Molecular Surface',
        params: MolecularSurfaceParams,
        get renderObjects() {
            const renderObjects = []
            if (currentProps.visuals.includes('surface')) renderObjects.push(...gaussianSurfaceRepr.renderObjects)
            if (currentProps.visuals.includes('wireframe')) renderObjects.push(...gaussianWireframeRepr.renderObjects)
            return renderObjects
        },
        get props() {
            return { ...gaussianSurfaceRepr.props, ...gaussianWireframeRepr.props, visuals: currentProps.visuals }
        },
        createOrUpdate: (props: Partial<MolecularSurfaceProps> = {}, structure?: Structure) => {
            if (structure) currentStructure = structure
            const qualityProps = getQualityProps(Object.assign({}, currentProps, props), currentStructure)
            currentProps = Object.assign({}, DefaultMolecularSurfaceProps, currentProps, props, qualityProps)
            return Task.create('Creating MolecularSurfaceRepresentation', async ctx => {
                if (currentProps.visuals.includes('surface')) await gaussianSurfaceRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
                if (currentProps.visuals.includes('wireframe')) await gaussianWireframeRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
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