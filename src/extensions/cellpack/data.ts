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
    location: 'surface' | 'interior' | 'cytoplasme'
    ingredients: Packing['ingredients']
    compartment?: CellCompartment
}

export interface CellCompartment {
    filename?: string
    geom_type?: 'raw' | 'file' | 'sphere' | 'mb' | 'None'
    compartment_primitives?: CompartmentPrimitives
}

export interface Cell {
    recipe: Recipe
    options?: RecipeOptions
    cytoplasme?: Packing
    compartments?: { [key: string]: Compartment }
    mapping_ids?: { [key: number]: [number, string] }
}

export interface RecipeOptions {
    resultfile?: string
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
    geom?: unknown
    geom_type?: 'raw' | 'file' | 'sphere' | 'mb' | 'None'
    mb?: CompartmentPrimitives
}

// Primitives discribing a compartment
export const enum CompartmentPrimitiveType {
    MetaBall = 0,
    Sphere = 1,
    Cube = 2,
    Cylinder = 3,
    Cone = 4,
    Plane = 5,
    None = 6
}

export interface CompartmentPrimitives{
    positions?: number[];
    radii?: number[];
    types?: CompartmentPrimitiveType[];
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
    principalVector?: Vec3;
    /** offset along membrane */
    offset?: Vec3;
    ingtype?: string;
    color?: Vec3;
    confidence?: number;
    Type?: string;
}

export interface IngredientSource {
    pdb: string;
    bu?: string; /** biological unit e.g AU,BU1,etc.. */
    selection?: string; /** NGL selection or :A or :B etc.. */
    model?: string; /** model number e.g 0,1,2... */
    transform: {
        center: boolean;
        translate?: Vec3;
    };
    biomt?: boolean;
}