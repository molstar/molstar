/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { Representation, RepresentationProps } from '..';
import { ColorTheme, SizeTheme } from '../../theme';
import { DefaultBaseProps, DefaultMeshProps } from '../util';

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export const DefaultStructureProps = {
    ...DefaultBaseProps,
    colorTheme: { name: 'instance-index' } as ColorTheme,
    sizeTheme: { name: 'physical' } as SizeTheme,
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