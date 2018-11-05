/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GaussianSurfaceVisual, GaussianSurfaceParams } from '../visual/gaussian-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { GaussianWireframeVisual, GaussianWireframeParams } from '../visual/gaussian-surface-wireframe';
import { paramDefaultValues, SelectParam } from 'mol-util/parameter';
import { GaussianDensityVolumeParams, GaussianDensityVolumeVisual } from '../visual/gaussian-density-volume';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const MolecularSurfaceParams = {
    ...GaussianSurfaceParams,
    ...GaussianWireframeParams,
    ...GaussianDensityVolumeParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
}
export const DefaultMolecularSurfaceProps = { ...paramDefaultValues(MolecularSurfaceParams), visuals: [ '0' ] }
export type MolecularSurfaceProps = typeof DefaultMolecularSurfaceProps

export type MolecularSurfaceRepresentation = StructureRepresentation<MolecularSurfaceProps>

export function MolecularSurfaceRepresentation(): MolecularSurfaceRepresentation {
    return Representation.createMulti('Molecular Surface', MolecularSurfaceParams, DefaultMolecularSurfaceProps, [
        UnitsRepresentation('Gaussian surface', GaussianSurfaceVisual),
        UnitsRepresentation('Gaussian wireframe', GaussianWireframeVisual),
        UnitsRepresentation('Gaussian volume', GaussianDensityVolumeVisual)
    ] as unknown as StructureRepresentation<MolecularSurfaceProps>[]) // TODO avoid cast to unknown
}