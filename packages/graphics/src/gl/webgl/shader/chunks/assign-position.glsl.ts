/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const assign_position = `
#ifdef dGeometryType_image
    mat4 transform = aTransform;
#else
    mat4 transform = applyTumble(aTransform);
#endif
#ifdef dGeometryType_textureMesh
    vec3 position = readFromTexture(tPosition, vertexId, uGeoTexDim).xyz;
#else
    vec3 position = aPosition;
#endif
#ifndef dGeometryType_image
    position = applyWiggle(position, float(group), int(aInstance));
#endif
#ifdef dGeometryType_mesh
    if (uScreenSize > 0.0) {
        // constant on screen: scaled about the instance's origin by the world size of a CSS pixel there, as
        // \`OverlayCamera\` bounds it
        float originDepth = max(-(uModelView * (uObjectTransform * transform[3])).z, 0.0);
        float pixelSize = (uIsOrtho == 1.0 ? 2.0 : 2.0 * originDepth) / (uProjection[1][1] * uCssHeight);
        position *= uScreenSize * pixelSize / uModelScale;
    }
#endif
// the position transformed, not the matrices multiplied first (three 4x4 products a vertex): the instance's
// transform, the object's transient one, then the model and view
vec4 placed = transform * vec4(position, 1.0);
vec4 world = uObjectTransform * placed;
// for accessing tColorGrid in vert shader and for clipping in frag shader
vModelPosition = (uModel * world).xyz;
// \`uModel\` is the camera's uniform scale (its scene group's view the identity), divided out as before
gridModelPosition = placed.xyz;
vec4 mvPosition = uModelView * world;
vViewPosition = mvPosition.xyz;
gl_Position = uProjection * mvPosition;
#ifdef dGeometryType_mesh
    if (uDepthOffset > 0.0) {
        // in depth only: its pixels stay; kept in front of the camera and within the depth range (as points'
        // \`uPointDepthOffset\`: closer than the offset, at the near plane)
        float offsetZ = mvPosition.z + uDepthOffset * uModelScale;
        if (uIsOrtho != 1.0) offsetZ = min(offsetZ, mvPosition.z * 1e-3);
        vec4 offsetPosition = uProjection * vec4(mvPosition.xy, offsetZ, 1.0);
        gl_Position.z = clamp(offsetPosition.z / offsetPosition.w, -1.0 + 1e-6, 1.0 - 1e-6) * gl_Position.w;
    }
#endif
`;
