/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { Representation, RepresentationProps } from '..';
import { ColorThemeProps, SizeThemeProps } from '../../theme';
import { DefaultBaseProps, DefaultMeshProps } from '../util';

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export const DefaultStructureProps = {
    ...DefaultBaseProps,
    colorTheme: { name: 'unit-index' } as ColorThemeProps,
    sizeTheme: { name: 'physical' } as SizeThemeProps,
}
export type StructureProps = typeof DefaultStructureProps

export const DefaultStructureMeshProps = {
    ...DefaultStructureProps,
    ...DefaultMeshProps
}
export type StructureMeshProps = typeof DefaultStructureMeshProps

export { ComplexRepresentation } from './complex-representation'
export { UnitsRepresentation } from './units-representation'
export { ComplexVisual } from './complex-visual'
export { UnitsVisual } from './units-visual'