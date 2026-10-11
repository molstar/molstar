/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { RenderableSchema, UniformBlock } from '../renderable/schema.js';
import {
    getElementKind,
    getStd140Members,
    getUniformBlock,
    getUniformBlockLayout,
    UniformBlockInfo,
    type UniformBlockLayout,
    type UniformBlockMember,
    UniformBlocks,
} from '../uniform-layout.js';
import type { UniformKind } from '../uniform.js';

/** Instance names of the blocks in WGSL, members are accessed as e.g. `uGlobal.uView` */
export const WgslUniformBlockInstance: { [k in UniformBlock]: string } = {
    global: 'uGlobal',
    material: 'uMaterial',
    object: 'uObject',
};

function getWgslElementType(kind: UniformKind): string {
    switch (kind) {
        case 'f':
            return 'f32';
        case 'i':
            return 'i32';
        // bool is not host-shareable, non-zero is true
        case 'b':
            return 'u32';
        case 'v2':
            return 'vec2<f32>';
        case 'v3':
            return 'vec3<f32>';
        case 'v4':
        case 'q':
            return 'vec4<f32>';
        case 'iv2':
            return 'vec2<i32>';
        case 'iv3':
            return 'vec3<i32>';
        case 'iv4':
            return 'vec4<i32>';
        case 'm3':
            return 'mat3x3<f32>';
        case 'm4':
            return 'mat4x4<f32>';
        default:
            throw new Error(`unsupported uniform block member kind '${kind}'`);
    }
}

/**
 * WGSL type matching the std140 layout of a member; array elements of 2-or-less
 * components are padded to 4-component vectors (uniform arrays need a 16-byte stride),
 * e.g. `f[]` becomes `array<vec4<f32>, N>` and element `i` is read as `a[i].x`.
 */
function getWgslType(m: UniformBlockMember): string {
    if (!m.length) return getWgslElementType(m.kind);
    const kind = getElementKind(m.kind);
    switch (kind) {
        case 'f':
        case 'v2':
            return `array<vec4<f32>, ${m.length}>`;
        case 'i':
        case 'iv2':
            return `array<vec4<i32>, ${m.length}>`;
        case 'b':
            return `array<vec4<u32>, ${m.length}>`;
        default:
            return `array<${getWgslElementType(kind)}, ${m.length}>`;
    }
}

/** WGSL struct and `var<uniform>` binding for `layout`, the member offsets equal the std140 layout */
export function getWgslUniformBlockDeclaration(layout: UniformBlockLayout, binding = 0) {
    const { name, bindingIndex } = UniformBlockInfo[layout.block];
    const lines = [`struct ${name} {`];
    for (const m of layout.members) {
        lines.push(`    ${m.name}: ${getWgslType(m)},`);
    }
    lines.push('};');
    lines.push(
        `@group(${bindingIndex}) @binding(${binding}) var<uniform> ${WgslUniformBlockInstance[layout.block]}: ${name};`,
    );
    return lines.join('\n');
}

/** Bind group of the per-draw uniforms (`Draw` struct, bound with a dynamic offset) */
export const WgslDrawGroup = 3;

/** std140 layout of the 'draw' uniforms of `schema`, written per draw */
export function getWgslDrawLayout(schema: RenderableSchema) {
    return getStd140Members('draw', schema);
}

/** WGSL struct and `var<uniform>` binding of the 'draw' uniforms, members are accessed as e.g. `uDraw.uDrawId` */
export function getWgslDrawDeclaration(schema: RenderableSchema) {
    const layout = getWgslDrawLayout(schema);
    if (!layout) return '';
    const lines = ['struct Draw {'];
    for (const m of layout.members) {
        lines.push(`    ${m.name}: ${getWgslType(m)},`);
    }
    lines.push('};');
    lines.push(`@group(${WgslDrawGroup}) @binding(0) var<uniform> uDraw: Draw;`);
    return lines.join('\n');
}

/**
 * WGSL declarations of the uniform blocks and the 'draw' uniforms of `schema`. WebGPU has no
 * classic uniforms, so the other uniforms without a block (arrays without `length`) are returned
 * as `unsupported`.
 */
export function getWgslUniformDeclarations(schema: RenderableSchema) {
    const lines: string[] = [];
    for (const block of UniformBlocks) {
        const layout = getUniformBlockLayout(block, schema);
        if (layout) lines.push(getWgslUniformBlockDeclaration(layout));
    }
    const draw = getWgslDrawDeclaration(schema);
    if (draw) lines.push(draw);
    const unsupported: string[] = [];
    for (const k of Object.keys(schema)) {
        const spec = schema[k];
        if (spec.type !== 'uniform' || spec.kind === 't' || spec.block === 'draw') continue;
        if (getUniformBlock(spec) === undefined) unsupported.push(k);
    }
    return { code: lines.join('\n\n'), unsupported };
}
