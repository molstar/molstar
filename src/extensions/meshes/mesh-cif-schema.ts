import * as MS from './molstar-lib-imports';

import Schema = MS.Column.Schema


const int = Schema.int;
const float = Schema.float;


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
export interface Mesh_Data_Database extends MS.Database<Mesh_Data_Schema> {}


// TODO in future, move to molstar/src/mol-io/reader/cif.ts: CIF.schema.mesh
export const CIF_schema_mesh = (frame: MS.CifFrame) => MS.toDatabase<Mesh_Data_Schema, Mesh_Data_Database>(Mesh_Data_Schema, frame);
