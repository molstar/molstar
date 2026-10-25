/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Items drawn as a quad of two triangles (6 vertices, no elements) with the item data pulled
 * from `tPositionGroup`: sets `mapping` (the quad corner, -1 or 1 per axis), `vertexId` (the item),
 * `position` and `group`
 */
export const assign_quad_position_group = `
vec2 mapping = vec2(1.0, 1.0); // vertices 2 and 5
#if __VERSION__ == 100
    int m = imod(VertexID, 6);
#else
    int m = VertexID % 6;
#endif
if (m == 0) {
    mapping = vec2(-1.0, 1.0);
} else if (m == 1 || m == 3) {
    mapping = vec2(-1.0, -1.0);
} else if (m == 4) {
    mapping = vec2(1.0, -1.0);
}

int vertexId = VertexID / 6;

vec4 positionGroup = read_tPositionGroup(vertexId);
vec3 position = positionGroup.rgb;
float group = positionGroup.a;
`;
