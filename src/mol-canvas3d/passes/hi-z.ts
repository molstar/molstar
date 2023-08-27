/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { HiZRenderable, createHiZRenderable } from '../../mol-gl/compute/hi-z';
import { isWebGL2 } from '../../mol-gl/webgl/compat';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Framebuffer } from '../../mol-gl/webgl/framebuffer';
import { Texture } from '../../mol-gl/webgl/texture';
import { fasterLog2 as _fasterLog2 } from '../../mol-math/approx';
import { Sphere3D } from '../../mol-math/geometry';
import { Mat4, Vec4 } from '../../mol-math/linear-algebra';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { arrayMinMax } from '../../mol-util/array';
import { isTimingMode } from '../../mol-util/debug';
import { ValueCell } from '../../mol-util/value-cell';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { DrawPass } from './draw';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const v4set = Vec4.set;
const fasterLog2 = _fasterLog2;

const MinLevel = 3;

function depthToViewZ(depth: number, near: number, far: number) {
    return (near * far) / ((far - near) * depth - far);
}

/**
 * Bounding rectangle of a clipped, perspective-projected 3D Sphere.
 * Michael Mara, Morgan McGuire. 2013
 *
 * Specialization by Arseny Kapoulkine, MIT License Copyright (c) 2018
 * https://github.com/zeux/niagara
 */
function projectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
    const prx = p[0] * r;
    const pry = p[1] * r;
    const prz = p[2] * r;
    const pzr2 = p[2] * p[2] - r * r;

    const vx = Math.sqrt(p[0] * p[0] + pzr2);
    const minx = ((vx * p[0] - prz) / (vx * p[2] + prx)) * projection[0];
    const maxx = ((vx * p[0] + prz) / (vx * p[2] - prx)) * projection[0];

    const vy = Math.sqrt(p[1] * p[1] + pzr2);
    const miny = ((vy * p[1] - prz) / (vy * p[2] + pry)) * projection[5];
    const maxy = ((vy * p[1] + prz) / (vy * p[2] - pry)) * projection[5];

    return v4set(
        out,
        1 - (maxx * 0.5 + 0.5),
        1 - (miny * -0.5 + 0.5),
        1 - (minx * 0.5 + 0.5),
        1 - (maxy * -0.5 + 0.5)
    );
}

type LevelData = {
    texture: Texture,
    framebuffer: Framebuffer,
    size: Vec2,
    invSize: Vec2,
    offset: number
}[]

export const HiZParams = {
    enabled: PD.Boolean(false, { description: 'Hierarchical Z-buffer occlusion culling. Only available for WebGL2.' }),
};
export type HiZProps = PD.Values<typeof HiZParams>

export class HiZPass {
    private enabled = false;

    private readonly camera = new Camera();
    private readonly aabb = Vec4();
    private readonly vp = Vec3();

    private readonly levelData: LevelData = [];
    private readonly fb: Framebuffer;
    private readonly buf: WebGLBuffer;
    private readonly tex: Texture;
    private readonly renderable: HiZRenderable;

    private sync: WebGLSync | null = null;
    private buffer = new Float32Array(0);
    private frame = 0;

    render(camera: Camera, props: HiZProps) {
        if (!props.enabled) return;

        const { gl, state } = this.webgl;
        if (!isWebGL2(gl) || this.sync !== null) return;

        this.camera.setState(camera.getSnapshot(), 0);
        Viewport.copy(this.camera.viewport, camera.viewport);

        if (isTimingMode) this.webgl.timer.mark('hi-Z');

        state.disable(gl.CULL_FACE);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.enable(gl.SCISSOR_TEST);
        state.depthMask(false);
        state.colorMask(true, true, true, true);
        state.clearColor(0, 0, 0, 0);

        this.fb.bind();
        this.tex.attachFramebuffer(this.fb, 0);

        const hw = this.tex.getWidth();
        const hh = this.tex.getHeight();

        state.viewport(0, 0, hw, hh);
        gl.clear(gl.COLOR_BUFFER_BIT);

        //

        const v = this.renderable.values;

        for (let i = 0, il = this.levelData.length; i < il; ++i) {
            const td = this.levelData[i];
            td.framebuffer.bind();

            ValueCell.update(v.uInvSize, td.invSize);
            if (i > 0) {
                ValueCell.update(v.tPreviousLevel, this.levelData[i - 1].texture);
            } else {
                ValueCell.update(v.tPreviousLevel, this.drawPass.depthTextureOpaque);
            }

            state.currentRenderItemId = -1;
            state.viewport(0, 0, td.size[0], td.size[1]);
            state.scissor(0, 0, td.size[0], td.size[1]);
            gl.clear(gl.COLOR_BUFFER_BIT);
            this.renderable.update();
            this.renderable.render();

            if (i >= MinLevel) {
                this.tex.bind(0);
                gl.copyTexSubImage2D(gl.TEXTURE_2D, 0, td.offset, 0, 0, 0, td.size[0], td.size[1]);
                this.tex.unbind(0);
            }
        }

        //

        this.fb.bind();
        this.tex.attachFramebuffer(this.fb, 0);

        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, this.buf);
        gl.bufferData(gl.PIXEL_PACK_BUFFER, this.buffer.byteLength, gl.STREAM_READ);
        gl.readPixels(camera.viewport.x, camera.viewport.y, hw, hh, gl.RED, gl.FLOAT, 0);
        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);

        this.sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
        gl.flush();

        if (isTimingMode) this.webgl.timer.markEnd('hi-Z');
    }

    tick(props: HiZProps) {
        const { gl } = this.webgl;
        if (!isWebGL2(gl)) return;

        this.enabled = props.enabled;
        if (!this.enabled || this.sync === null) return;

        const res = gl.clientWaitSync(this.sync, 0, 0);
        if (res === gl.WAIT_FAILED) {
            // console.log(`failed to get buffer data after ${this.frame + 1} frames`);
            gl.deleteSync(this.sync);
            this.sync = null;
            this.frame = 0;
        } else if (res === gl.TIMEOUT_EXPIRED) {
            this.frame += 1;
            // console.log(`waiting for buffer data for ${this.frame} frames`);
        } else {
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, this.buf);
            gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, this.buffer);
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
            // console.log(`got buffer data after ${this.frame + 1} frames`);
            gl.deleteSync(this.sync);
            this.sync = null;
            this.frame = 0;

            // const p = PixelData.flipY(PixelData.create(this.dest.slice(), this.tex.getWidth(), this.tex.getHeight()));
            // for (let i = 0, il = p.array.length; i < il; ++i) {
            //     p.array[i] = -p.array[i];
            // }
            // console.log('fff', arrayMinMax(this.dest), arrayMinMax(p.array));
            // printTextureImage(p, { scale: 3, id: 'foo', useCanvas: true, normalize: true, pixelated: true });

            this.camera.update();
        }
    }

    isOccluded = (s: Sphere3D) => {
        if (!this.enabled) return false;

        const { camera, vp, aabb } = this;
        const { near, far, viewport, view, projection } = camera;

        const r = (s.radius * 1.2) + 1.52;
        v3transformMat4(vp, s.center, view);
        projectSphere(aabb, vp, r, projection);
        const w = aabb[2] - aabb[0];
        const h = aabb[3] - aabb[1];

        const pr = Math.max(w * viewport.width, h * viewport.height);
        const lod = Math.ceil(fasterLog2(pr / 2));
        if (lod >= this.levelData.length || lod < MinLevel) return false;

        const { offset, size } = this.levelData[lod];

        const u = aabb[0] + w / 2;
        const v = aabb[1] + h / 2;

        const ts = size[0];
        const x = u * ts;
        const y = v * ts;
        const z = vp[2] + r;

        const dx = Math.floor(x);
        const dy = Math.ceil(y);
        const dw = this.tex.getWidth();

        if (dx + 1 >= ts || dy + 1 >= ts) return false;

        const di = (ts - dy - 1) * dw + dx + offset;
        if (z > depthToViewZ(this.buffer[di], near, far)) return false;

        const di1 = (ts - dy - 1) * dw + dx + 1 + offset;
        if (z > depthToViewZ(this.buffer[di1], near, far)) return false;

        const di2 = (ts - dy + 1 - 1) * dw + dx + offset;
        if (z > depthToViewZ(this.buffer[di2], near, far)) return false;

        const di3 = (ts - dy + 1 - 1) * dw + dx + 1 + offset;
        if (z > depthToViewZ(this.buffer[di3], near, far)) return false;

        return true;
    };

    syncSize() {
        if (!this.webgl.isWebGL2) return;

        const width = this.drawPass.colorTarget.getWidth();
        const height = this.drawPass.colorTarget.getHeight();

        const levels = Math.ceil(Math.log(Math.max(width, height)) / Math.log(2));
        this.buffer = new Float32Array(Math.pow(2, levels - MinLevel) * Math.pow(2, levels - 1 - MinLevel));
        this.tex.define(Math.pow(2, levels - MinLevel), Math.pow(2, levels - 1 - MinLevel));
        this.levelData.length = 0;
        for (let i = 0; i < levels; ++i) {
            const framebuffer = this.webgl.resources.framebuffer();
            const levelSize = Math.pow(2, levels - i - 1);
            const size = Vec2.create(levelSize, levelSize);
            const invSize = Vec2.create(1 / levelSize, 1 / levelSize);
            const texture = this.webgl.resources.texture('image-float32', 'alpha', 'float', 'nearest');
            texture.define(levelSize, levelSize);
            texture.attachFramebuffer(framebuffer, 0);
            this.levelData.push({ texture, framebuffer, size, invSize, offset: 0 });
        }

        let offset = 0;
        for (let i = 0, il = levels; i < il; ++i) {
            const td = this.levelData[i];
            if (i >= MinLevel) {
                this.levelData[i].offset = offset;
                offset += td.size[0];
            }
        }
    }

    getDebugContext() {
        return {
            levelData: this.levelData,
            buffer: this.buffer,
            tex: this.tex,
            camera: this.camera,
        };
    }

    constructor(private webgl: WebGLContext, private drawPass: DrawPass) {
        if (!webgl.isWebGL2) {
            this.enabled = false;
            return;
        }

        const buf = webgl.gl.createBuffer();
        if (buf === null) throw new Error('Could not create WebGL buffer');

        this.fb = webgl.resources.framebuffer();
        this.buf = buf;
        this.tex = webgl.resources.texture('image-float32', 'alpha', 'float', 'nearest');
        this.renderable = createHiZRenderable(webgl, this.drawPass.depthTextureOpaque);
    }
}

//

export class HiZDebug {
    private readonly container = document.createElement('div');
    private readonly canvas = document.createElement('canvas');

    private aabb = Vec4();
    private vp = Vec3();

    getRect(element: Element) {
        return document.getElementById('molstar-hiz-rect') || this.addRect(element);
    }

    addRect(element: Element) {
        const rect = document.createElement('div');
        rect.id = 'molstar-hiz-rect';

        Object.assign(rect.style, {
            display: 'none',
            pointerEvents: 'none',
        });

        element.parentElement?.appendChild(rect);

        return rect;
    }

    updateRect(p: Vec4, o: boolean, element: Element, viewport: Viewport) {
        const { x, y, width, height } = viewport;
        const rect = this.getRect(element);

        const minx = p[0] * width + x;
        const miny = p[1] * height + y;
        const maxx = p[2] * width + x;
        const maxy = p[3] * height + y;

        Object.assign(rect.style, {
            border: o ? 'solid red' : 'solid green',
            display: 'block',
            position: 'absolute',
            left: `${minx}px`,
            top: `${miny}px`,
            width: `${maxx - minx}px`,
            height: `${maxy - miny}px`,
        });
    }

    printOcclusion(s: Sphere3D, element: Element) {
        const { camera, levelData, tex, buffer } = this.hiz.getDebugContext();
        const { vp, aabb } = this;
        const { near, far, viewport, view, projection } = camera;

        const r = (s.radius * 1.2) + 1.52;
        v3transformMat4(vp, s.center, view);
        projectSphere(aabb, vp, r, projection);
        const w = aabb[2] - aabb[0];
        const h = aabb[3] - aabb[1];

        const pr = Math.max(w * viewport.width, h * viewport.height);
        const lod = Math.ceil(fasterLog2(pr / 2));
        if (lod >= levelData.length || lod < MinLevel) return false;
        const td = levelData[lod];
        const { offset } = td;

        vp[2] *= -1;
        vp[2] -= r;

        const u = aabb[0] + w / 2;
        const v = aabb[1] + h / 2;

        const x = u * td.size[0];
        const y = v * td.size[1];

        const dx = Math.floor(x);
        const dy = Math.floor(y);
        const dw = tex.getWidth();
        const dh = tex.getHeight();
        const ds = dw * dh;
        const tw = td.size[0];
        const th = td.size[1];

        let occluded = true;
        if (dx + 1 >= tw || dy + 1 >= th) occluded = false;

        const di = (th - dy - 1) * dw + dx + offset;
        if (vp[2] <= -depthToViewZ(buffer[di], near, far)) occluded = false;

        const dz = buffer[di];
        const dd = -(near * far) / ((far - near) * buffer[di] - far);
        console.log({ depth: dd, objectDepth: vp[2], z: dz }, arrayMinMax(buffer), { near, far, lod, tw, th }, { offset, dx, dy, di, dw, dh, ds }, { x, y, u, v, w, h }, aabb, td);
        const data = new Uint8ClampedArray(tw * th * 4);
        data.fill(255);
        let min = +Infinity;
        let max = -Infinity;
        let mind = +Infinity;
        let maxd = -Infinity;
        for (let y = 0; y < th; ++y) {
            for (let x = 0; x < tw; ++x) {
                const i = (th - y - 1) * (tw) + x;
                const dv = buffer[y * dw + x + offset];
                const v = dv * 255;
                data[i * 4 + 0] = v;
                data[i * 4 + 3] = (255 - v) / 2;
                min = Math.min(min, v);
                max = Math.max(max, v);
                mind = Math.min(mind, dv);
                maxd = Math.max(maxd, dv);
            }
        }
        const imageData = new ImageData(data, tw, th);
        console.log({ min, max, mind, maxd, tw, th });
        this.print(imageData, element);
        this.updateRect(aabb, occluded, element, camera.viewport);

        return occluded;
    }

    private print(imageData: ImageData, element: Element) {
        Object.assign(this.container.style, {
            position: 'absolute',
            bottom: '0px',
            top: '0px',
            left: '0px',
            right: '0px',
            pointerEvents: 'none',
        });
        element.parentElement?.appendChild(this.container);

        this.canvas.width = imageData.width;
        this.canvas.height = imageData.height;
        const outCtx = this.canvas.getContext('2d');
        if (!outCtx) throw new Error('Could not create canvas 2d context');
        outCtx.putImageData(imageData, 0, 0);
        this.canvas.style.width = '100%';
        this.canvas.style.height = '100%';
        this.canvas.style.imageRendering = 'pixelated';
        this.canvas.style.position = 'relative';
        this.canvas.style.pointerEvents = 'none';
        // this.canvas.style.opacity = '0.5';
        this.container.appendChild(this.canvas);
    }

    constructor(private readonly hiz: HiZPass) {

    }
}
