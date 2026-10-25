/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Vertices per page of the vertex arena of chunked isosurfaces (a triangle may span two pages) */
export const IsosurfaceChunkSize = 4096;

/**
 * Words of a page table before its page ids: vertex count, vertex page count, status, drawn count, the min (3 f32) and
 * the step (3 f32) of the box the surface's positions are quantized in, triangle count, triangle page count, vertex
 * capacity (see `mc_chunk_common_wgsl`)
 */
export const IsosurfaceTableHeader = 13;

/** Bytes per vertex of the arena, and per triangle (3 indices): 3 words, see `iso_vertex_wgsl` */
export const IsosurfaceVertexBytes = 12;
