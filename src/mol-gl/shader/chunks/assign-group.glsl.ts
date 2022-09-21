

import { loadCheckpoint } from '../../../mol-util/debug';
loadCheckpoint(`mol-gl/shader/chunks/assign-group.glsl.ts::start`);
export const assign_group = `
#ifdef dGeometryType_textureMesh
    float group = unpackRGBToInt(readFromTexture(tGroup, VertexID, uGeoTexDim).rgb);
#else
    float group = aGroup;
#endif
`;
loadCheckpoint(`mol-gl/shader/chunks/assign-group.glsl.ts::end`);
