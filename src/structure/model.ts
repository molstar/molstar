/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data'
import { Vec3, Mat4 } from '../utils/linear-algebra'


/** Unit = essentially a list of residues (usually a chain) */
export interface Unit extends Readonly<{
    /** The static part (essentially residue annotations) */
    structre: Unit.Structure,
    /** 3D arrangement that often changes with time. */
    conformation: Unit.Conformation
}> { }

export namespace Unit {
    export interface Structure extends Readonly<{
        data: Data.Structure,
        /** A globally unique number for this instance (to easily determine unique structures within a model) */
        key: number,
        /** Reference to the data.entities table */
        entity: number,
        /** Reference to the data.chains table */
        chain: number,
        /** Indices into the data.residues table. */
        residues: ArrayLike<number>,
        /** Offsets of atoms in the residue layer. start = offsets[i], endExclusive = offsets[i + 1] */
        atomOffsets: ArrayLike<number>,
        /** Index of a residue in the corresponding residues array. */
        atomResidue: number
    }> { }

    export interface Bonds extends Readonly<{
        /**
         * Where bonds for atom A start and end.
         * Start at idx, end at idx + 1
         */
        offset: ArrayLike<number>,
        neighbor: ArrayLike<number>,

        order: ArrayLike<number>,
        flags: ArrayLike<number>,

        count: number
    }> { }

    export interface Conformation extends Readonly<{
        positions: Data.Positions
        spatialLookup: any, // TODO
        boundingSphere: { readonly center: Vec3, readonly radius: number },
        secodaryStructure: Data.SecondaryStructures,
        bonds: Bonds
    }> { }

    export type OperatorKind =
        | { kind: 'identity' }
        | { kind: 'symmetry', id: string, hkl: Vec3 }
        | { kind: 'assembly', index: number }

    export interface Operator extends Readonly<{
        kind: OperatorKind,
        transform: Mat4,
        inverse: Mat4
    }> { }

    export interface Lookup3D {
        // TODO
    }
}

export interface Model extends Readonly<{
    index: number,

    structure: Model.Structure,
    conformation: Model.Conformation
}> { }

export namespace Model {
    // TODO: model residues between units
    // use a map of map of bonds?

    export interface Structure extends Readonly<{
        operators: Unit.Operator[],
        units: Unit.Structure[]
    }> { }

    export interface Conformation extends Readonly<{
        units: Unit.Conformation[],
        spatialLookup: Unit.Lookup3D
    }> { }
}

export namespace Atom {
    /**
     * Represents a "packed reference" to an atom.
     * This is because selections can then be represented
     * a both number[] and Float64Array(), making it much
     * more efficient than storing an array of objects.
     */
    export type Reference = number
    export interface Location { unit: number, atom: number }

    const { _uint32, _float64 } = (function() {
        const data = new ArrayBuffer(8);
        return { _uint32: new Uint32Array(data), _float64: new Float64Array(data) };
    }());

    export function emptyLocation(): Location { return { unit: 0, atom: 0 }; }

    export function getRef(location: Location) {
        _uint32[0] = location.unit;
        _uint32[1] = location.atom;
        return _float64[0];
    }

    export function getLocation(ref: Reference) {
        return updateLocation(ref, emptyLocation());
    }

    export function updateLocation(ref: Reference, location: Location): Location {
        _float64[0] = ref;
        location.unit = _uint32[0];
        location.atom = _uint32[1];
        return location;
    }
}
