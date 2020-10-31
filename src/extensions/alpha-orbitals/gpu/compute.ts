/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { QuadSchema, QuadValues } from '../../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../../mol-gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, Values } from '../../../mol-gl/renderable/schema';
import { ShaderCode } from '../../../mol-gl/shader-code';
import quad_vert from '../../../mol-gl/shader/quad.vert';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../../mol-gl/webgl/render-item';
import { ValueCell } from '../../../mol-util';
import { arrayMin } from '../../../mol-util/array';
import { isLittleEndian } from '../../../mol-util/is-little-endian';
import { CollocationParams } from '../collocation';
import { normalizeBasicOrder } from '../orbitals';
import shader_frag from './shader.frag';

const AlphaOrbitalsSchema = {
    ...QuadSchema,
    uDimensions: UniformSpec('v3'),
    uMin: UniformSpec('v3'),
    uDelta: UniformSpec('v3'),
    tCenters: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    tInfo: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    tCoeff: TextureSpec('image-float32', 'rgb', 'float', 'nearest'),
    tAlpha: TextureSpec('image-float32', 'alpha', 'float', 'nearest'),
    uNCenters: UniformSpec('i'),
    uNAlpha: UniformSpec('i'),
    uNCoeff: UniformSpec('i'),
    uMaxCoeffs: UniformSpec('i'),
    uLittleEndian: UniformSpec('i') // TODO: boolean uniforms
};
const AlphaOrbitalsShaderCode = ShaderCode('postprocessing', quad_vert, shader_frag);
type AlphaOrbitalsRenderable = ComputeRenderable<Values<typeof AlphaOrbitalsSchema>>

function createTextureData({
    basis,
    sphericalOrder,
    alphaOrbitals,
    cutoffThreshold
}: CollocationParams) {
    let centerCount = 0;
    let baseCount = 0;
    let coeffCount = 0;
    for (const atom of basis.atoms) {
        for (const shell of atom.shells) {
            for (const L of shell.angularMomentum) {
                if (L > 4) {
                    // TODO: will L > 4 be required? Would need to precompute more functions in that case.
                    throw new Error('Angular momentum L > 4 not supported.');
                }

                centerCount++;
                baseCount += 2 * L + 1;
                coeffCount += shell.exponents.length;
            }
        }
    }

    const centers = new Float32Array(4 * centerCount);
    // L, alpha_offset, coeff_offset_start, coeff_offset_end
    const info = new Float32Array(4 * centerCount);
    const alpha = new Float32Array(baseCount);
    const coeff = new Float32Array(3 * coeffCount);

    let maxCoeffs = 0;
    let cO = 0, aO = 0, coeffO = 0;
    for (const atom of basis.atoms) {
        for (const shell of atom.shells) {

            let amIndex = 0;
            for (const L of shell.angularMomentum) {
                const a0 = normalizeBasicOrder(L, alphaOrbitals.slice(aO, aO + 2 * L + 1), sphericalOrder);

                const cutoffRadius = cutoffThreshold > 0
                    ? Math.sqrt(-Math.log(cutoffThreshold) / arrayMin(shell.exponents))
                    : 10000;

                centers[4 * cO + 0] = atom.center[0];
                centers[4 * cO + 1] = atom.center[1];
                centers[4 * cO + 2] = atom.center[2];
                centers[4 * cO + 3] = cutoffRadius * cutoffRadius;

                info[4 * cO + 0] = L;
                info[4 * cO + 1] = aO;
                info[4 * cO + 2] = coeffO;
                info[4 * cO + 3] = coeffO + shell.exponents.length;

                for (let i = 0; i < a0.length; i++) alpha[aO + i] = a0[i];

                const c0 = shell.coefficients[amIndex++];
                for (let i = 0; i < shell.exponents.length; i++) {
                    coeff[3 * (coeffO + i) + 0] = c0[i];
                    coeff[3 * (coeffO + i) + 1] = shell.exponents[i];
                }

                if (c0.length > maxCoeffs) {
                    maxCoeffs = c0.length;
                }

                cO++;
                aO += 2 * L + 1;
                coeffO += shell.exponents.length;
            }
        }
    }

    return { nCenters: centerCount, nAlpha: baseCount, nCoeff: coeffCount, maxCoeffs, centers, info, alpha, coeff };
}

function getPostprocessingRenderable(ctx: WebGLContext, params: CollocationParams): AlphaOrbitalsRenderable {
    const data = createTextureData(params);

    const values: Values<typeof AlphaOrbitalsSchema> = {
        ...QuadValues,
        uDimensions: ValueCell.create(params.grid.dimensions),
        uMin: ValueCell.create(params.grid.box.min),
        uDelta: ValueCell.create(params.grid.delta),
        uNCenters: ValueCell.create(data.nCenters),
        uNAlpha: ValueCell.create(data.nAlpha),
        uNCoeff: ValueCell.create(data.nCoeff),
        uMaxCoeffs: ValueCell.create(data.maxCoeffs),
        tCenters: ValueCell.create({ width: data.nCenters, height: 1, array: data.centers }),
        tInfo: ValueCell.create({ width: data.nCenters, height: 1, array: data.info }),
        tCoeff: ValueCell.create({ width: data.nCoeff, height: 1, array: data.coeff }),
        tAlpha: ValueCell.create({ width: data.nAlpha, height: 1, array: data.alpha }),
        uLittleEndian: ValueCell.create(isLittleEndian()),
    };

    const schema = { ...AlphaOrbitalsSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', AlphaOrbitalsShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function normalizeParams(webgl: WebGLContext) {
    if (!webgl.isWebGL2) {
        // workaround for webgl1 limitation that loop counters need to be `const`
        (AlphaOrbitalsSchema.uNCenters as any) = DefineSpec('number');
        (AlphaOrbitalsSchema.uMaxCoeffs as any) = DefineSpec('number');
    }
}

export function gpuComputeAlphaOrbitalsGridValues(webgl: WebGLContext, params: CollocationParams) {
    const [nx, ny, nz] = params.grid.dimensions;

    normalizeParams(webgl);

    if (!webgl.computeTargets['alpha-oribtals']) {
        webgl.computeTargets['alpha-oribtals'] = webgl.createRenderTarget(nx, ny * nz, false, 'uint8', 'nearest');
    } else {
        webgl.computeTargets['alpha-oribtals'].setSize(nx, ny * nz);
    }

    const target = webgl.computeTargets['alpha-oribtals'];
    const renderable = getPostprocessingRenderable(webgl, params);

    const width = nx;
    const height = ny * nz;

    const { gl, state } = webgl;
    target.bind();
    gl.viewport(0, 0, width, height);
    gl.scissor(0, 0, width, height);
    state.disable(gl.SCISSOR_TEST);
    state.disable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.depthMask(false);
    renderable.render();

    const array = new Uint8Array(width * height * 4);
    webgl.readPixels(0, 0, width, height, array);
    const floats = new Float32Array(array.buffer, array.byteOffset, width * height);
    renderable.dispose();

    return floats;
}

export function canComputeAlphaOrbitalsOnGPU(webgl?: WebGLContext) {
    return !!webgl?.extensions.textureFloat;
}