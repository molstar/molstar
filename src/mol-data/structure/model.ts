/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './data'
import { Selectors } from './selectors'
import { Vec3, Mat4 } from '../../mol-base/math/linear-algebra'

let _uid = 0;
/** Model-related unique identifiers */
export function getUid() { return _uid++; }

/** Unit = essentially a list of residues (usually a chain) */
export interface Unit extends Readonly<{
    operator: Unit.Operator,
    /** The static part (essentially residue annotations) */
    structure: Unit.Structure,
    /** 3D arrangement that often changes with time. */
    conformation: Unit.Conformation,
    /** Position of i-th atom. Special function for this because it's the only one that depends on "operator" */
    atomPosition(i: number, p: Vec3): Vec3,
    /** Property access */
    selectors: Selectors
}> { }

export namespace Unit {
    export interface Structure extends Readonly<{
        /** A globally unique identifier for this instance (to easily determine unique structures within a model) */
        id: number,
        /** Source data for this structure */
        data: Data.Structure,
        /** Reference to the data.entities table */
        entity: number,
        /** Reference to the data.chains table */
        chain: number,
        /** Indices into the data.residues table. */
        residues: ArrayLike<number>,
        /** Offsets of atoms in the residue layer. start = offsets[i], endExclusive = offsets[i + 1] */
        atomOffsets: ArrayLike<number>,
        /** Indices into data.atoms table */
        atoms: ArrayLike<number>
        /** Index of a residue in the corresponding residues array. */
        atomResidue: ArrayLike<number>
    }> { }

    export interface Conformation extends Readonly<{
        /** A globally unique identifier for this instance (to easily determine unique structures within a model) */
        id: number,

        positions: Data.Positions,

        /** Per residue secondary structure assignment. */
        secondaryStructure: Data.SecondaryStructures
    }> {
        '@spatialLookup': any, // TODO
        '@boundingSphere': Readonly<{ center: Vec3, radius: number }>,
        '@bonds': Data.Bonds
    }

    export namespace Conformation {
        export function spatialLookup(conformation: Conformation): any {
            throw 'not implemented'
        }

        export function boundingSphere(conformation: Conformation): any {
            throw 'not implemented'
        }

        export function bonds(conformation: Conformation): any {
            throw 'not implemented'
        }
    }

    export type OperatorKind =
        | { kind: 'identity' }
        | { kind: 'symmetry', id: string, hkl: Vec3 }
        | { kind: 'assembly', index: number }
        | { kind: 'custom' }

    export interface Operator extends Readonly<{
        kind: OperatorKind,
        transform: Mat4,
        inverse: Mat4
    }> { }

    export interface SpatialLookup {
        // TODO
    }
}

export interface Model extends Readonly<{
    units: Unit[],

    structure: { [id: number]: Unit.Structure },
    conformation: { [id: number]: Unit.Conformation }
}> {
    '@unitLookup'?: Unit.SpatialLookup,
    /** bonds between units */
    '@unitBonds'?: Data.Bonds
}

export namespace Model {
    export function unitLookup(model: Model): Unit.SpatialLookup {
        throw 'not implemented';
    }

    export function unitBonds(model: Model): Data.Bonds {
        throw 'not implemented';
    }

    export function join(a: Model, b: Model): Model {
        return {
            units: [...a.units, ...b.units],
            structure: { ...a.structure, ...b.structure },
            conformation: { ...a.conformation, ...b.conformation }
        }
    }
}

export namespace Atom {
    /**
     * Represents a "packed reference" to an atom.
     * This is because selections can then be represented
     * a both number[] and Float64Array(), making it much
     * more efficient than storing an array of objects.
     */
    export type PackedReference = number
    export interface Reference { unit: number, atom: number }
    export type Selection = ArrayLike<PackedReference>

    const { _uint32, _float64 } = (function() {
        const data = new ArrayBuffer(8);
        return { _uint32: new Uint32Array(data), _float64: new Float64Array(data) };
    }());

    export function emptyRef(): Reference { return { unit: 0, atom: 0 }; }

    export function packRef(location: Reference) {
        _uint32[0] = location.unit;
        _uint32[1] = location.atom;
        return _float64[0];
    }

    export function unpackRef(ref: PackedReference) {
        return updateRef(ref, emptyRef());
    }

    export function updateRef(ref: PackedReference, location: Reference): Reference {
        _float64[0] = ref;
        location.unit = _uint32[0];
        location.atom = _uint32[1];
        return location;
    }
}

export interface Selector<T> {
    '@type': T,
    create(m: Model, u: number): Selector.Field<T>
}

export function Selector<T>(create: (m: Model, u: number) => (a: number, v?: T) => T): Selector<T> {
    return { create } as Selector<T>;
}

export namespace Selector {
    export type Field<T> = (a: number, v?: T) => T;

    export type Category = { [s: string]: Selector<any> };
    export type Unit<C extends Category> = { [F in keyof C]: Field<C[F]['@type']> }

    export type Set<S extends { [n: string]: Category }> = { [C in keyof S]: Unit<S[C]> }
}