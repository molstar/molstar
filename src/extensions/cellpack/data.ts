/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Quat } from '../../mol-math/linear-algebra';

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
    paths: [string, string][] // [name: string, path: string][]
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

export interface Positions {
    coords?: Vec3[];
}

export interface Radii {
    radii?: number[];
}

export interface Ingredient {
    source: IngredientSource;
    results: [Vec3, Quat][];
    name: string;
    /** Vec3[]];CoarseGraind Beads coordinates LOD */
    positions?: [Positions];
    /** number[]];CoarseGraind Beads radii LOD */
    radii?: [Radii];
    /** Number of `curveX` properties in the object where `X` is a 0-indexed number */
    nbCurve?: number;
    uLength?: number;
    /** Curve properties are Vec3[] but that is not expressable in TypeScript */
    [curveX: string]: unknown;
    /** the orientation in the membrane */
    principalAxis?: Vec3;
    /** offset along membrane */
    offset?: Vec3;
    ingtype?: string;
    color?: Vec3;
    confidence?: number;
}

export interface IngredientSource {
    pdb: string;
    bu?: string;  /** biological unit e.g AU,BU1,etc.. */
    selection?: string; /** NGL selection or :A or :B etc.. */
    model?: string;     /** model number e.g 0,1,2... */
    transform: {
        center: boolean;
        translate?: Vec3;
    };
    biomt?: boolean;
}