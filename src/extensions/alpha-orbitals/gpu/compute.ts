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
import { Basis } from '../cubes';
import { normalizeBasicOrder, SphericalBasisOrder } from '../orbitals';
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
    uWidth: UniformSpec('f'),
    uNCenters: UniformSpec('i'),
    uNAlpha: UniformSpec('i'),
    uNCoeff: UniformSpec('i'),
    uMaxCoeffs: UniformSpec('i'),
    uLittleEndian: UniformSpec('b'),

    uOccupancy: UniformSpec('f'),
    tCumulativeSum: TextureSpec('texture', 'rgba', 'ubyte', 'nearest')
};
type AlphaOrbitalsSchema = Values<typeof AlphaOrbitalsSchema>
const AlphaOrbitalsName = 'alpha-orbitals';
const AlphaOrbitalsShaderCode = ShaderCode(AlphaOrbitalsName, quad_vert, shader_frag);
type AlphaOrbitalsRenderable = ComputeRenderable<AlphaOrbitalsSchema>

function getNormalizedAlpha(basis: Basis, alphaOrbitals: number[], sphericalOrder: SphericalBasisOrder) {
    const alpha = new Float32Array(alphaOrbitals.length);

    let aO = 0;
    for (const atom of basis.atoms) {
        for (const shell of atom.shells) {
            for (const L of shell.angularMomentum) {
                const a0 = normalizeBasicOrder(L, alphaOrbitals.slice(aO, aO + 2 * L + 1), sphericalOrder);
                for (let i = 0; i < a0.length; i++) alpha[aO + i] = a0[i];
                aO += 2 * L + 1;
            }
        }
    }

    return alpha;
}

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

function createAlphaOrbitalsRenderable(ctx: WebGLContext, params: CollocationParams): AlphaOrbitalsRenderable {
    const data = createTextureData(params);

    const [nx, ny, nz] = params.grid.dimensions;
    const width = Math.ceil(Math.sqrt(nx * ny * nz));

    if (!ctx.namedTextures[AlphaOrbitalsName]) {
        ctx.namedTextures[AlphaOrbitalsName] = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }

    const values: AlphaOrbitalsSchema = {
        ...QuadValues,
        uDimensions: ValueCell.create(params.grid.dimensions),
        uMin: ValueCell.create(params.grid.box.min),
        uDelta: ValueCell.create(params.grid.delta),
        uWidth: ValueCell.create(width),
        uNCenters: ValueCell.create(data.nCenters),
        uNAlpha: ValueCell.create(data.nAlpha),
        uNCoeff: ValueCell.create(data.nCoeff),
        uMaxCoeffs: ValueCell.create(data.maxCoeffs),
        tCenters: ValueCell.create({ width: data.nCenters, height: 1, array: data.centers }),
        tInfo: ValueCell.create({ width: data.nCenters, height: 1, array: data.info }),
        tCoeff: ValueCell.create({ width: data.nCoeff, height: 1, array: data.coeff }),
        tAlpha: ValueCell.create({ width: data.nAlpha, height: 1, array: data.alpha }),
        uLittleEndian: ValueCell.create(isLittleEndian()),

        uOccupancy: ValueCell.create(0),
        tCumulativeSum: ValueCell.create(ctx.namedTextures[AlphaOrbitalsName])
    };

    const schema = { ...AlphaOrbitalsSchema };
    if (!ctx.isWebGL2) {
        // workaround for webgl1 limitation that loop counters need to be `const`
        (schema.uNCenters as any) = DefineSpec('number');
        (schema.uMaxCoeffs as any) = DefineSpec('number');
    }

    const renderItem = createComputeRenderItem(ctx, 'triangles', AlphaOrbitalsShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function getAlphaOrbitalsRenderable(ctx: WebGLContext, params: CollocationParams): AlphaOrbitalsRenderable {
    if (ctx.namedComputeRenderables[AlphaOrbitalsName]) {
        const v = ctx.namedComputeRenderables[AlphaOrbitalsName].values as AlphaOrbitalsSchema;

        const data = createTextureData(params);

        const [nx, ny, nz] = params.grid.dimensions;
        const width = Math.ceil(Math.sqrt(nx * ny * nz));

        ValueCell.update(v.uDimensions, params.grid.dimensions);
        ValueCell.update(v.uMin, params.grid.box.min);
        ValueCell.update(v.uDelta, params.grid.delta);
        ValueCell.updateIfChanged(v.uWidth, width);
        ValueCell.updateIfChanged(v.uNCenters, data.nCenters);
        ValueCell.updateIfChanged(v.uNAlpha, data.nAlpha);
        ValueCell.updateIfChanged(v.uNCoeff, data.nCoeff);
        ValueCell.updateIfChanged(v.uMaxCoeffs, data.maxCoeffs);
        ValueCell.update(v.tCenters, { width: data.nCenters, height: 1, array: data.centers });
        ValueCell.update(v.tInfo, { width: data.nCenters, height: 1, array: data.info });
        ValueCell.update(v.tCoeff, { width: data.nCoeff, height: 1, array: data.coeff });
        ValueCell.update(v.tAlpha, { width: data.nAlpha, height: 1, array: data.alpha });
        ValueCell.updateIfChanged(v.uLittleEndian, isLittleEndian());
        ValueCell.updateIfChanged(v.uOccupancy, 0);

        ctx.namedComputeRenderables[AlphaOrbitalsName].update();
    } else {
        ctx.namedComputeRenderables[AlphaOrbitalsName] = createAlphaOrbitalsRenderable(ctx, params);
    }
    return ctx.namedComputeRenderables[AlphaOrbitalsName];
}

export function gpuComputeAlphaOrbitalsGridValues(webgl: WebGLContext, params: CollocationParams) {
    const [nx, ny, nz] = params.grid.dimensions;
    const renderable = getAlphaOrbitalsRenderable(webgl, params);
    const width = renderable.values.uWidth.ref.value;

    if (!webgl.namedFramebuffers[AlphaOrbitalsName]) {
        webgl.namedFramebuffers[AlphaOrbitalsName] = webgl.resources.framebuffer();
    }
    const framebuffer = webgl.namedFramebuffers[AlphaOrbitalsName];

    if (!webgl.namedTextures[AlphaOrbitalsName]) {
        webgl.namedTextures[AlphaOrbitalsName] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }
    webgl.namedTextures[AlphaOrbitalsName].define(width, width);
    webgl.namedTextures[AlphaOrbitalsName].attachFramebuffer(framebuffer, 'color0');

    const { gl, state } = webgl;
    framebuffer.bind();
    gl.viewport(0, 0, width, width);
    gl.scissor(0, 0, width, width);
    state.disable(gl.SCISSOR_TEST);
    state.disable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.depthMask(false);
    renderable.render();

    const array = new Uint8Array(width * width * 4);
    webgl.readPixels(0, 0, width, width, array);
    return new Float32Array(array.buffer, array.byteOffset, nx * ny * nz);
}

export function canComputeAlphaOrbitalsOnGPU(webgl?: WebGLContext) {
    return !!webgl?.extensions.textureFloat;
}

const AlphaOrbitalsDensity0 = AlphaOrbitalsName;
const AlphaOrbitalsDensity1 = AlphaOrbitalsName + '1';

export function gpuComputeAlphaOrbitalsDensityGridValues(webgl: WebGLContext, params: CollocationParams, orbitals: number[][]) {
    const [nx, ny, nz] = params.grid.dimensions;
    const renderable = getAlphaOrbitalsRenderable(webgl, params);
    const width = renderable.values.uWidth.ref.value;

    if (!webgl.namedFramebuffers[AlphaOrbitalsDensity0]) {
        webgl.namedFramebuffers[AlphaOrbitalsDensity0] = webgl.resources.framebuffer();
    }
    const framebuffer = webgl.namedFramebuffers[AlphaOrbitalsDensity0];

    if (!webgl.namedTextures[AlphaOrbitalsDensity0]) {
        webgl.namedTextures[AlphaOrbitalsDensity0] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }
    if (!webgl.namedTextures[AlphaOrbitalsDensity1]) {
        webgl.namedTextures[AlphaOrbitalsDensity1] = webgl.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }


    const tex = [webgl.namedTextures[AlphaOrbitalsDensity0], webgl.namedTextures[AlphaOrbitalsDensity1]];

    tex[0].define(width, width);
    tex[1].define(width, width);

    const values = renderable.values as AlphaOrbitalsSchema;
    const { gl, state } = webgl;

    gl.viewport(0, 0, width, width);
    gl.scissor(0, 0, width, width);
    state.disable(gl.SCISSOR_TEST);
    state.disable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.depthMask(false);

    gl.clearColor(0, 0, 0, 0);

    tex[0].attachFramebuffer(framebuffer, 'color0');
    gl.clear(gl.COLOR_BUFFER_BIT);

    tex[1].attachFramebuffer(framebuffer, 'color0');
    gl.clear(gl.COLOR_BUFFER_BIT);

    for (let i = 0; i < 108; i++) {
        const alpha = getNormalizedAlpha(params.basis, orbitals[i], params.sphericalOrder);

        ValueCell.update(values.uOccupancy, 2);
        ValueCell.update(values.tCumulativeSum, tex[(i + 1) % 2]);
        ValueCell.update(values.tAlpha, { width: alpha.length, height: 1, array: alpha });
        renderable.update();
        tex[i % 2].attachFramebuffer(framebuffer, 'color0');
        renderable.render();
    }

    const array = new Uint8Array(width * width * 4);
    webgl.readPixels(0, 0, width, width, array);

    tex[0].define(1, 1);
    tex[1].define(1, 1);

    return new Float32Array(array.buffer, array.byteOffset, nx * ny * nz);
}