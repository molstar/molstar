/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Quat } from '../../../../mol-math/linear-algebra';

export interface CellPack {
    cell: Cell
    packings: CellPacking[]
}

export interface CellPacking {
    name: string,
    location: 'surface' | 'interior' | 'cytoplasme',
    ingredients: Packing['ingredients']
}

//

export interface Cell {
    recipe: Recipe
    cytoplasme?: Packing
    compartments?: { [key: string]: Compartment }
}

export interface Recipe {
    setupfile: string
    /** First entry is name, secound is path: [name: string, path: string][] */
    paths: [string, string][]
    version: string
    name: string
}

export interface Compartment {
    surface?: Packing
    interior?: Packing
}

export interface Packing {
    ingredients: { [key: string]: Ingredient }
}

export interface Ingredient {
    source: IngredientSource
    results: [Vec3, Quat][]
    name: string
    positions?: [Vec3[]] // why wrapped in an extra array?
    radii?: [number[]] // why wrapped in an extra array?

    /** Number of `curveX` properties in the object where `X` is a 0-indexed number */
    nbCurve?: number
    /** Curve properties are Vec3[] but that is not expressable in TypeScript */
    [curveX: string]: unknown
}

export interface IngredientSource {
    pdb: string
    transform: { center: boolean, translate?: Vec3 }
    biomt?: boolean
}
