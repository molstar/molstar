/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Column, Database } from '../../mol-data/db';
import { CifFrame } from '../../mol-io/reader/cif';
import { toDatabase } from '../../mol-io/reader/cif/schema';


const int = Column.Schema.int;
const float = Column.Schema.float;


// TODO in future, move to molstar/src/mol-io/reader/cif/schema/mesh.ts
export const Mesh_Data_Schema = {
    mesh: {
        id: int,
    },
    mesh_vertex: {
        mesh_id: int,
        vertex_id: int,
        x: float,
        y: float,
        z: float,
    },
    /** Table of triangles, 3 rows per triangle */
    mesh_triangle: {
        mesh_id: int,
        /** Indices of vertices within mesh */
        vertex_id: int,
    }
};
export type Mesh_Data_Schema = typeof Mesh_Data_Schema;
export interface Mesh_Data_Database extends Database<Mesh_Data_Schema> {}


// TODO in future, move to molstar/src/mol-io/reader/cif.ts: CIF.schema.mesh
export const CIF_schema_mesh = (frame: CifFrame) => toDatabase<Mesh_Data_Schema, Mesh_Data_Database>(Mesh_Data_Schema, frame);
