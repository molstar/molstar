/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MoleculeType, ComponentType } from '../types'

export interface ChemicalComponent {
    id: string
    type: ComponentType
    moleculeType: MoleculeType
    name: string
    synonyms: string[]
    formula: string
    formulaWeight: number
}

export type ChemicalComponentMap = ReadonlyMap<string, ChemicalComponent>