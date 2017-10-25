/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO: define property accessor intefaces, graphs, spatial lookups and what have you.

interface Model {

    // Incremented when data changes
    dataVersion: number,

    // Incremented when the underlying conformation changes
    conformationVersion: number,
}

export default Model
