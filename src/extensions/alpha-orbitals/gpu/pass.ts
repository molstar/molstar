/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { QuadSchema, QuadValues } from '../../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../../mol-gl/renderable';
import { TextureSpec, UniformSpec, Values } from '../../../mol-gl/renderable/schema';
import { ShaderCode } from '../../../mol-gl/shader-code';
import quad_vert from '../../../mol-gl/shader/quad.vert';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../../mol-gl/webgl/render-item';
import { RenderTarget } from '../../../mol-gl/webgl/render-target';
import { ValueCell } from '../../../mol-util';
import { CollocationParams } from '../collocation';
import { normalizeBasicOrder } from '../orbitals';
import shader_frag from './shader.frag';

const AlphaOrbitalsSchema = {
    ...QuadSchema,
    uDimensions: UniformSpec('v3'),
    uMin: UniformSpec('v3'),
    uDelta: UniformSpec('v3'),
    tCenters: TextureSpec('image-float32', 'rgb', 'float', 'nearest'),
    tInfo: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    tCoeff: TextureSpec('image-float32', 'rgb', 'float', 'nearest'),
    tAlpha: TextureSpec('image-float32', 'alpha', 'float', 'nearest'),
    uNCenters: UniformSpec('i'),
    uNAlpha: UniformSpec('i'),
    uNCoeff: UniformSpec('i'),
    uLittleEndian: UniformSpec('i') // TODO: boolean uniforms
};
const AlphaOrbitalsShaderCode = ShaderCode('postprocessing', quad_vert, shader_frag);
type AlphaOrbitalsRenderable = ComputeRenderable<Values<typeof AlphaOrbitalsSchema>>

function createTextureData({
    basis,
    sphericalOrder,
    alphaOrbitals,
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

    const centers = new Float32Array(3 * centerCount);
    // L, alpha_offset, coeff_offset_start, coeff_offset_end
    const info = new Float32Array(4 * centerCount);
    const alpha = new Float32Array(baseCount);
    const coeff = new Float32Array(3 * coeffCount);

    let cO = 0, aO = 0, coeffO = 0;
    for (const atom of basis.atoms) {
        for (const shell of atom.shells) {

            let amIndex = 0;
            for (const L of shell.angularMomentum) {
                const a0 = normalizeBasicOrder(L, alphaOrbitals.slice(aO, aO + 2 * L + 1), sphericalOrder);

                if (cO === 1) {
                    console.log('y', atom.center[1]);
                }

                centers[3 * cO + 0] = atom.center[0];
                centers[3 * cO + 1] = atom.center[1];
                centers[3 * cO + 2] = atom.center[2];

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

                cO++;
                aO += 2 * L + 1;
                coeffO += shell.exponents.length;
            }
        }
    }

    return { nCenters: centerCount, nAlpha: baseCount, nCoeff: coeffCount, centers, info, alpha, coeff };
}

function getPostprocessingRenderable(ctx: WebGLContext, params: CollocationParams): AlphaOrbitalsRenderable {
    const data = createTextureData(params);

    console.log(data);

    const values: Values<typeof AlphaOrbitalsSchema> = {
        ...QuadValues,
        uDimensions: ValueCell.create(params.grid.dimensions),
        uMin: ValueCell.create(params.grid.box.min),
        uDelta: ValueCell.create(params.grid.delta),
        uNCenters: ValueCell.create(data.nCenters),
        uNAlpha: ValueCell.create(data.nAlpha),
        uNCoeff: ValueCell.create(data.nCoeff),
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

export class AlphaOrbitalsPass {
    target: RenderTarget
    renderable: AlphaOrbitalsRenderable

    constructor(private webgl: WebGLContext, private params: CollocationParams) {
        const [nx, ny, nz] = params.grid.dimensions;

        // TODO: add single component float32 render target option for WebGL2?
        // TODO: figure out the ordering so that it does not have to be remapped in the shader
        this.target = webgl.createRenderTarget(nx, ny * nz, false, 'uint8', 'nearest');
        this.renderable = getPostprocessingRenderable(webgl, params);
    }

    render() {
        const [nx, ny, nz] = this.params.grid.dimensions;
        const width = nx;
        const height = ny * nz;

        // const { x, y, width, height } = this.camera.viewport;

        const { gl, state } = this.webgl;
        this.target.bind();
        // this.setQuadShift(0, 0);
        gl.viewport(0, 0, width, height);
        gl.scissor(0, 0, width, height);
        state.disable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
        this.renderable.render();
    }

    getData() {
        const [nx, ny, nz] = this.params.grid.dimensions;
        const width = nx;
        const height = ny * nz;

        this.render();
        this.target.bind();
        const array = new Uint8Array(width * height * 4);
        this.webgl.readPixels(0, 0, width, height, array);
        // PixelData.flipY({ array, width, height });
        const floats = new Float32Array(array.buffer, array.byteOffset, width * height);

        // console.log(array);
        // console.log(floats);

        return floats;

        // return new ImageData(new Uint8ClampedArray(array), width, height);
    }
}

function isLittleEndian() {
    const arrayBuffer = new ArrayBuffer(2);
    const uint8Array = new Uint8Array(arrayBuffer);
    const uint16array = new Uint16Array(arrayBuffer);
    uint8Array[0] = 0xAA; // set first byte
    uint8Array[1] = 0xBB; // set second byte
    if(uint16array[0] === 0xBBAA) return 1;
    return 0;
}