/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RenderableSchema, Values, KindValue, UniformSpec, TextureSpec, DefineSpec, RenderableValues } from '../renderable/schema';
import { WebGLContext } from '../webgl/context';
import { getRegularGrid3dDelta, RegularGrid3d } from '../../mol-math/geometry/common';
import shader_template from '../shader/util/grid3d-template.frag';
import quad_vert from '../shader/quad.vert';
import { ShaderCode } from '../shader-code';
import { UUID, ValueCell } from '../../mol-util';
import { objectForEach } from '../../mol-util/object';
import { getUniformGlslType, isUniformValueScalar } from '../webgl/uniform';
import { QuadSchema, QuadValues } from './util';
import { createComputeRenderItem } from '../webgl/render-item';
import { createComputeRenderable } from '../renderable';
import { isLittleEndian } from '../../mol-util/is-little-endian';
import { RuntimeContext } from '../../mol-task';

export interface Grid3DComputeRenderableSpec<S extends RenderableSchema, P, CS> {
    schema: S,
    // indicate which params are loop bounds for WebGL1 compat
    loopBounds?: (keyof S)[]
    utilCode?: string,
    mainCode: string,
    returnCode: string,

    values(params: P, grid: RegularGrid3d): SchemaValues<S>,

    cumulative?: {
        init(params: P, values: SchemaValues<S>): CS,
        next(params: P, state: CS, values: Values<S>): boolean
    }
}

type SchemaValues<S extends RenderableSchema> = { readonly [k in keyof S]: KindValue[S[k]['kind']] }

const FrameBufferName = 'grid3d-computable' as const;
const Texture0Name = 'grid3d-computable-0' as const;
const Texture1Name = 'grid3d-computable-1' as const;

const SchemaBase = {
    ...QuadSchema,
    uDimensions: UniformSpec('v3'),
    uMin: UniformSpec('v3'),
    uDelta: UniformSpec('v3'),
    uWidth: UniformSpec('f'),
    uLittleEndian: UniformSpec('b'),
};

const CumulativeSumSchema = {
    tCumulativeSum: TextureSpec('texture', 'rgba', 'ubyte', 'nearest')
};

export function createGrid3dComputable<S extends RenderableSchema, P, CS>(spec: Grid3DComputeRenderableSpec<S, P, CS>) {
    const id = UUID.create22();

    const uniforms: string[] = [];

    objectForEach(spec.schema, (u, k) => {
        if (u.type === 'define') return;
        if (u.kind.indexOf('[]') >= 0) throw new Error('array uniforms are not supported');
        const isBound = (spec.loopBounds?.indexOf(k) ?? -1) >= 0;
        if (isBound) uniforms.push(`#ifndef ${k}`);
        uniforms.push(`uniform ${getUniformGlslType(u.kind as any)} ${k};`);
        if (isBound) uniforms.push(`#endif`);
    });

    const code = shader_template
        .replace('{UNIFORMS}', uniforms.join('\n'))
        .replace('{UTILS}', spec.utilCode ?? '')
        .replace('{MAIN}', spec.mainCode)
        .replace('{RETURN}', spec.returnCode);

    const shader = ShaderCode(id, quad_vert, code);

    return async (ctx: RuntimeContext, webgl: WebGLContext, grid: RegularGrid3d, params: P) => {
        const schema: RenderableSchema = {
            ...SchemaBase,
            ...(spec.cumulative ? CumulativeSumSchema : {}),
            ...spec.schema,
        };

        if (!webgl.isWebGL2) {
            if (spec.loopBounds) {
                for (const b of spec.loopBounds) {
                    (schema as any)[b] = DefineSpec('number');
                }
            }
            (schema as any)['WEBGL1'] = DefineSpec('boolean');
        }

        if (spec.cumulative) {
            (schema as any)['CUMULATIVE'] = DefineSpec('boolean');
        }

        if (!webgl.namedFramebuffers[FrameBufferName]) {
            webgl.namedFramebuffers[FrameBufferName] = webgl.resources.framebuffer();
        }
        const framebuffer = webgl.namedFramebuffers[FrameBufferName];

        if (!webgl.namedTextures[Texture0Name]) {
            webgl.namedTextures[Texture0Name] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        }
        if (spec.cumulative && !webgl.namedTextures[Texture1Name]) {
            webgl.namedTextures[Texture1Name] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        }

        const tex = [webgl.namedTextures[Texture0Name], webgl.namedTextures[Texture1Name]];

        const [nx, ny, nz] = grid.dimensions;
        const uWidth = Math.ceil(Math.sqrt(nx * ny * nz));

        const values: SchemaValues<S & typeof SchemaBase> = {
            ...QuadValues,
            uDimensions: grid.dimensions,
            uMin: grid.box.min,
            uDelta: getRegularGrid3dDelta(grid),
            uWidth,
            uLittleEndian: isLittleEndian(),
            ...spec.values(params, grid)
        } as any;

        if (!webgl.isWebGL2) {
            (values as any).WEBGL1 = true;
        }
        if (spec.cumulative) {
            (values as any).tCumulativeSum = tex[0];
            (values as any).CUMULATIVE = true;
        }

        let renderable = webgl.namedComputeRenderables[id];
        if (renderable) {
            const cells = renderable.values as RenderableValues;
            objectForEach(values, (c, k) => {
                const s = schema[k];
                if (s?.type === 'value' || s?.type === 'attribute') return;

                if (!s || !isUniformValueScalar(s.kind as any)) {
                    ValueCell.update(cells[k], c);
                } else {
                    ValueCell.updateIfChanged(cells[k], c);
                }
            });
        } else {
            const cells = {} as any;
            objectForEach(values, (v, k) => cells[k] = ValueCell.create(v));
            renderable = createComputeRenderable(createComputeRenderItem(webgl, 'triangles', shader, schema, cells), cells);
        }

        if (spec.cumulative) {
            throw new Error('nyi');
        } else {
            tex[0].attachFramebuffer(framebuffer, 'color0');
            framebuffer.bind();
            resetGl(webgl, uWidth);
            renderable.render();
        }

        const array = new Uint8Array(uWidth * uWidth * 4);
        webgl.readPixels(0, 0, uWidth, uWidth, array);
        return new Float32Array(array.buffer, array.byteOffset, nx * ny * nz);
    };
}

function resetGl(webgl: WebGLContext, w: number) {
    const { gl, state } = webgl;
    gl.viewport(0, 0, w, w);
    gl.scissor(0, 0, w, w);
    state.disable(gl.SCISSOR_TEST);
    state.disable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.depthMask(false);
}