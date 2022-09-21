

import { loadCheckpoint } from '../../../mol-util/debug';
loadCheckpoint(`mol-gl/shader/chunks/assign-position.glsl.ts::start`);
export const assign_position = `
mat4 model = uModel * aTransform;
mat4 modelView = uView * model;
#ifdef dGeometryType_textureMesh
    vec3 position = readFromTexture(tPosition, VertexID, uGeoTexDim).xyz;
#else
    vec3 position = aPosition;
#endif
vec4 position4 = vec4(position, 1.0);
// for accessing tColorGrid in vert shader and for clipping in frag shader
vModelPosition = (model * position4).xyz;
vec4 mvPosition = modelView * position4;
vViewPosition = mvPosition.xyz;
gl_Position = uProjection * mvPosition;
`;
loadCheckpoint(`mol-gl/shader/chunks/assign-position.glsl.ts::end`);
