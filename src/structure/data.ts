/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type Table<Data> = { rowCount: number } & { [E in keyof Data]: ArrayLike<Data[E]> }
export type DataTable<Data> = { data: any } & Table<Data & { dataIndex: number }>

export interface Position { x: number, y: number, z: number }
export interface Positions extends Table<Position> {}

export interface Atom {
    name: string,
    elementSymbol: string,
    altLoc: string | null,
}
export interface Atoms extends DataTable<Atom> { }

export interface Residue {
    key: number,
    name: string,
    seqNumber: number,
    insCode: string | null,
    isHet: number
}
export interface Residues extends DataTable<Residue> { }

export interface Chain { key: number, id: string }
export interface Chains extends DataTable<Chain> { }

export interface Entity { key: number, id: string }
export interface Entities extends DataTable<Entity> { }

export type SourceData =
    | { kind: 'mmCIF', data: any  } // TODO
    | { kind: 'custom', data: any  } // TODO

export interface Structure {
    atoms: Atoms,
    residues: Residues,
    chains: Chains,
    entities: Entities
}

export interface SecondaryStructure {

}
export interface SecondaryStructures extends Table<SecondaryStructure> {

}