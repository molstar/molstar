/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// import { File as GroFile } from 'mol-io/reader/gro/schema'
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif'
import CIF, { CifFrame } from 'mol-io/reader/cif';
import { PdbFile } from 'mol-io/reader/pdb/schema';

type Format =
    // | Format.gro
    | Format.mmCIF

namespace Format {
    // export interface gro { kind: 'gro', data: GroFile }
    export interface mmCIF { kind: 'mmCIF', data: mmCIF_Database, frame: CifFrame }
    export function mmCIF(frame: CifFrame, data?: mmCIF_Database): mmCIF { return { kind: 'mmCIF', data: data || CIF.schema.mmCIF(frame), frame }; }

    export interface PDB { kind: 'PDB', data: PdbFile }
    export function PDB(data: PdbFile) { return { kind: 'PDB', data }; }
}

export default Format