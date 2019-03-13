/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Cage {
    readonly vertices: ArrayLike<number>
    readonly edges: ArrayLike<number>
}

export function createCage(vertices: ArrayLike<number>, edges: ArrayLike<number>): Cage {
    return { vertices, edges }
}