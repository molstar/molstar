/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueIndex } from 'mol-model/structure';
import { BuiltInStructureRepresentationsName } from 'mol-repr/structure/registry';
import { BuiltInColorThemeName } from 'mol-theme/color';

export interface StructureInfo {
    ligands: { name: string, indices: ResidueIndex[] }[],
    assemblies: { id: string, description: string, isPreferred: boolean }[]
}

export type SupportedFormats = 'cif' | 'pdb'
export interface LoadParams {
    url: string,
    format?: SupportedFormats,
    assemblyId?: string,
    representationStyle?: RepresentationStyle
}

export interface RepresentationStyle {
    sequence?: RepresentationStyle.Entry,
    hetGroups?: RepresentationStyle.Entry,
    water?: RepresentationStyle.Entry
}

export namespace RepresentationStyle {
    export type Entry = { kind?: BuiltInStructureRepresentationsName, coloring?: BuiltInColorThemeName }
}