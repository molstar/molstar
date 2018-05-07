/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;

uniform float alpha;

#pragma glslify: import('./chunks/color-frag-params.glsl')

void main(){
    #pragma glslify: import('./chunks/color-assign-material.glsl')
    gl_FragColor = vec4(material, alpha);
}