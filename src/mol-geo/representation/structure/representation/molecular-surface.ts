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
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { PickingId } from '../../../geometry/picking';
import { Task } from 'mol-task';
import { GaussianWireframeVisual, GaussianWireframeParams } from '../visual/gaussian-surface-wireframe';
import { getQualityProps } from '../../util';
import { paramDefaultValues, MultiSelectParam, SelectParam } from 'mol-view/parameter';
import { GaussianDensityVolumeParams, GaussianDensityVolumeVisual } from '../visual/gaussian-density-volume';
import { SizeThemeName, SizeThemeOptions } from 'mol-view/theme/size';

const VisualOptions = [['surface', 'Surface'], ['wireframe', 'Wireframe'], ['volume', 'Volume']] as [string, string][]

export const MolecularSurfaceParams = {
    ...GaussianSurfaceParams,
    ...GaussianWireframeParams,
    ...GaussianDensityVolumeParams,

    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
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
    const gaussianVolumeRepr = UnitsRepresentation('Gaussian volume', GaussianDensityVolumeVisual)
    return {
        label: 'Molecular Surface',
        params: MolecularSurfaceParams,
        get renderObjects() {
            const renderObjects = []
            if (currentProps.visuals.includes('surface')) renderObjects.push(...gaussianSurfaceRepr.renderObjects)
            if (currentProps.visuals.includes('wireframe')) renderObjects.push(...gaussianWireframeRepr.renderObjects)
            if (currentProps.visuals.includes('volume')) renderObjects.push(...gaussianVolumeRepr.renderObjects)
            return renderObjects
        },
        get props() {
            return { ...gaussianSurfaceRepr.props, ...gaussianWireframeRepr.props, ...gaussianVolumeRepr.props, visuals: currentProps.visuals }
        },
        createOrUpdate: (props: Partial<MolecularSurfaceProps> = {}, structure?: Structure) => {
            if (structure) currentStructure = structure
            const qualityProps = getQualityProps(Object.assign({}, currentProps, props), currentStructure)
            currentProps = Object.assign({}, DefaultMolecularSurfaceProps, currentProps, props, qualityProps)
            return Task.create('Creating MolecularSurfaceRepresentation', async ctx => {
                if (currentProps.visuals.includes('surface')) await gaussianSurfaceRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
                if (currentProps.visuals.includes('wireframe')) await gaussianWireframeRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
                if (currentProps.visuals.includes('volume')) await gaussianVolumeRepr.createOrUpdate(currentProps, currentStructure).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const surfaceLoci = gaussianSurfaceRepr.getLoci(pickingId)
            const wireframeLoci = gaussianWireframeRepr.getLoci(pickingId)
            const volumeLoci = gaussianVolumeRepr.getLoci(pickingId)
            if (isEmptyLoci(surfaceLoci)) {
                if (isEmptyLoci(wireframeLoci)) {
                    return volumeLoci
                } else {
                    return wireframeLoci
                }
            } else {
                return surfaceLoci
            }
        },
        mark: (loci: Loci, action: MarkerAction) => {
            const markSurfaceElement = gaussianSurfaceRepr.mark(loci, action)
            const markWireframeElement = gaussianWireframeRepr.mark(loci, action)
            const markVolumeElement = gaussianVolumeRepr.mark(loci, action)
            return markSurfaceElement || markWireframeElement || markVolumeElement
        },
        destroy() {
            gaussianSurfaceRepr.destroy()
            gaussianWireframeRepr.destroy()
            gaussianVolumeRepr.destroy()
        }
    }
}