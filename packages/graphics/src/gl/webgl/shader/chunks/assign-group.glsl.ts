/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const assign_group = `
#if defined(dInstanceData)
    readInstanceData();
#endif
#ifdef dGeometryType_textureMesh
    GroupId group = unpackGroupId(readFromTexture(tGroup, vertexId, uGeoTexDim));
#else
    GroupId group = GroupId(aGroup);
#endif
`;
