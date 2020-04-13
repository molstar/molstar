/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export default `
float matrixScale(in mat4 m){
    vec4 r = m[0];
    return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}
`;