/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import Format from './format'
import StructureSequence from './properties/sequence'
import { AtomicHierarchy, AtomicConformation } from './properties/atomic'
import { ModelSymmetry } from './properties/symmetry'
import { CoarseHierarchy, CoarseConformation } from './properties/coarse'
import { Entities } from './properties/common';
import { SecondaryStructure } from './properties/seconday-structure';

// import from_gro from './formats/gro'
import from_mmCIF from './formats/mmcif'
/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: UUID,
    label: string,

    modelNum: number,

    sourceData: Format,

    symmetry: ModelSymmetry,
    entities: Entities,
    sequence: StructureSequence,

    atomicHierarchy: AtomicHierarchy,
    atomicConformation: AtomicConformation,

    /** Various parts of the code can "cache" custom properties here */
    properties: {
        readonly secondaryStructure: SecondaryStructure,
        // maps modified residue name to its parent
        readonly modifiedResidueNameMap: Map<string, string>,
        readonly asymIdSerialMap: Map<string, number>,
        [customName: string]: any
    },

    // TODO: separate properties to "static" (propagated with trajectory) and "dynamic" (computed for each frame separately)

    coarseHierarchy: CoarseHierarchy,
    coarseConformation: CoarseConformation
}> {

} { }

namespace Model {
    export function create(format: Format) {
        switch (format.kind) {
            // case 'gro': return from_gro(format);
            case 'mmCIF': return from_mmCIF(format);
        }
    }

    // TODO: figure the place to include this?
    // export interface Property<T, K> { (model: Model, index: number): T, _kind: K }
    // export interface AtomicProperty<T> extends Property<T, 'atomic'> { }
    // export interface CoarseProperty<T> extends Property<T, 'coarse'> { }
    // export interface SphereProperty<T> extends Property<T, 'sphere'> { }
    // export interface GaussianProperty<T> extends Property<T, 'gaussian'> { }

    // export function atomProp<T>(p: (model: Model, i: number) => T): AtomicProperty<T> { return p as any; }
    // export function residueProp<T>(p: (model: Model, residueIndex: number) => T): AtomicProperty<T> {
    //     return p as any;
    // }
}

export default Model