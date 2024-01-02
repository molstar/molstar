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
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { ValueCell } from '../../mol-util/value-cell';
import { Camera } from '../camera';
import { Viewport } from '../camera/util';
import { DrawPass } from './draw';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { getBuffer } from '../../mol-gl/webgl/buffer';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const v4set = Vec4.set;
const fasterLog2 = _fasterLog2;

function perspectiveDepthToViewZ(depth: number, near: number, far: number) {
    return (near * far) / ((far - near) * depth - far);
}

function orthographicDepthToViewZ(depth: number, near: number, far: number) {
    return depth * (near - far) - near;
}

function depthToViewZ(depth: number, near: number, far: number, projection: Mat4) {
    return projection[11] === -1
        ? perspectiveDepthToViewZ(depth, near, far)
        : orthographicDepthToViewZ(depth, near, far);
}

/**
 * Bounding rectangle of a clipped, perspective-projected 3D Sphere.
 * Michael Mara, Morgan McGuire. 2013
 *
 * Specialization by Arseny Kapoulkine, MIT License Copyright (c) 2018
 * https://github.com/zeux/niagara
 */
function perspectiveProjectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
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
        maxx * -0.5 + 0.5,
        miny * 0.5 + 0.5,
        minx * -0.5 + 0.5,
        maxy * 0.5 + 0.5
    );
}

function orthographicProjectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
    const sx = projection[0];
    const sy = projection[5];

    const minx = (p[0] + r) * sx;
    const maxx = (p[0] - r) * sx;
    const miny = (p[1] + r) * sy;
    const maxy = (p[1] - r) * sy;

    return v4set(
        out,
        maxx * 0.5 + 0.5,
        miny * -0.5 + 0.5,
        minx * 0.5 + 0.5,
        maxy * -0.5 + 0.5
    );
}

function projectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
    return projection[11] === -1
        ? perspectiveProjectSphere(out, p, r, projection)
        : orthographicProjectSphere(out, p, r, projection);
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
    maxFrameLag: PD.Numeric(10, { min: 1, max: 30, step: 1 }, { description: 'Maximum number of frames to wait for Z-buffer data.' }),
    minLevel: PD.Numeric(3, { min: 1, max: 10, step: 1 }),
};
export type HiZProps = PD.Values<typeof HiZParams>

export class HiZPass {
    private readonly viewport = Viewport();

    private near = 0;
    private far = 0;
    private readonly view = Mat4();
    private readonly projection = Mat4();

    private nextNear = 0;
    private nextFar = 0;
    private readonly nextView = Mat4();
    private readonly nextProjection = Mat4();

    private readonly aabb = Vec4();
    private readonly vp = Vec3();

    private readonly levelData: LevelData = [];
    private readonly fb: Framebuffer;
    private readonly buf: WebGLBuffer;
    private readonly tex: Texture;
    private readonly renderable: HiZRenderable;
    private readonly supported: boolean;

    private sync: WebGLSync | null = null;
    private buffer = new Float32Array(0);
    private frameLag = 0;
    private ready = false;

    readonly props: HiZProps;

    clear() {
        if (!this.supported) return;

        const { gl } = this.webgl;
        if (!isWebGL2(gl)) return;

        if (this.sync !== null) {
            gl.deleteSync(this.sync);
            this.sync = null;
        }
        this.frameLag = 0;
        this.ready = false;

        if (this.debug) {
            this.debug.rect.style.display = 'none';
            this.debug.container.style.display = 'none';
        }
    }

    render(camera: Camera) {
        if (!this.supported || !this.props.enabled) return;

        const { gl, state } = this.webgl;
        if (!isWebGL2(gl) || this.sync !== null) return;

        this.nextNear = camera.near;
        this.nextFar = camera.far;
        Mat4.copy(this.nextView, camera.view);
        Mat4.copy(this.nextProjection, camera.projection);

        if (isTimingMode) this.webgl.timer.mark('hi-Z');

        state.disable(gl.CULL_FACE);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.disable(gl.SCISSOR_TEST);
        state.depthMask(false);
        state.colorMask(true, true, true, true);
        state.clearColor(0, 0, 0, 0);

        //

        const v = this.renderable.values;
        const s = Math.pow(2, Math.ceil(Math.log(Math.max(gl.drawingBufferWidth, gl.drawingBufferHeight)) / Math.log(2)) - 1);

        for (let i = 0, il = this.levelData.length; i < il; ++i) {
            const td = this.levelData[i];
            td.framebuffer.bind();

            if (i > 0) {
                ValueCell.update(v.uInvSize, td.invSize);
                ValueCell.update(v.uOffset, Vec2.set(v.uOffset.ref.value, 0, 0));
                ValueCell.update(v.tPreviousLevel, this.levelData[i - 1].texture);
            } else {
                ValueCell.update(v.uInvSize, Vec2.set(v.uInvSize.ref.value, 1 / s, 1 / s));
                ValueCell.update(v.uOffset, Vec2.set(v.uOffset.ref.value,
                    this.viewport.x / gl.drawingBufferWidth,
                    this.viewport.y / gl.drawingBufferHeight
                ));
                ValueCell.update(v.tPreviousLevel, this.drawPass.depthTextureOpaque);
            }

            state.currentRenderItemId = -1;
            state.viewport(0, 0, td.size[0], td.size[1]);
            gl.clear(gl.COLOR_BUFFER_BIT);
            this.renderable.update();
            this.renderable.render();

            if (i >= this.props.minLevel) {
                this.tex.bind(0);
                gl.copyTexSubImage2D(gl.TEXTURE_2D, 0, td.offset, 0, 0, 0, td.size[0], td.size[1]);
                this.tex.unbind(0);
            }
        }

        //

        this.tex.attachFramebuffer(this.fb, 0);
        const hw = this.tex.getWidth();
        const hh = this.tex.getHeight();

        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, this.buf);
        gl.bufferData(gl.PIXEL_PACK_BUFFER, this.buffer.byteLength, gl.STREAM_READ);
        gl.readPixels(0, 0, hw, hh, gl.RED, gl.FLOAT, 0);
        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);

        this.sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
        gl.flush();

        if (isTimingMode) this.webgl.timer.markEnd('hi-Z');
    }

    tick() {
        if (!this.supported || !this.props.enabled || this.sync === null) return;

        const { gl } = this.webgl;
        if (!isWebGL2(gl)) return;

        const res = gl.clientWaitSync(this.sync, 0, 0);
        if (res === gl.WAIT_FAILED || this.frameLag >= this.props.maxFrameLag) {
            // console.log(`failed to get buffer data after ${this.frameLag + 1} frames`);
            gl.deleteSync(this.sync);
            this.sync = null;
            this.frameLag = 0;
            this.ready = false;
        } else if (res === gl.TIMEOUT_EXPIRED) {
            this.frameLag += 1;
            // console.log(`waiting for buffer data for ${this.frameLag} frames`);
        } else {
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, this.buf);
            gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, this.buffer);
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
            // console.log(`got buffer data after ${this.frameLag + 1} frames`);
            gl.deleteSync(this.sync);
            this.sync = null;
            this.frameLag = 0;
            this.ready = true;

            // if (isDebugMode) {
            //     const p = PixelData.flipY(PixelData.create(this.buffer.slice(), this.tex.getWidth(), this.tex.getHeight()));
            //     printTextureImage(p, { scale: MinLevel, id: 'hiz', useCanvas: true, pixelated: true });
            // }

            this.near = this.nextNear;
            this.far = this.nextFar;
            Mat4.copy(this.view, this.nextView);
            Mat4.copy(this.projection, this.nextProjection);
        }
    }

    private transform(s: Sphere3D) {
        const { view, vp } = this;
        v3transformMat4(vp, s.center, view);
        const r = (s.radius * 1.2) + 1.52;
        return { vp, r };
    }

    private project(vp: Vec3, r: number) {
        const { projection, aabb, viewport } = this;
        projectSphere(aabb, vp, r, projection);
        const w = aabb[2] - aabb[0];
        const h = aabb[3] - aabb[1];
        const pr = Math.max(w * viewport.width, h * viewport.height);
        const lod = Math.ceil(fasterLog2(pr / 2));
        return { aabb, w, h, pr, lod };
    }

    isOccluded = (s: Sphere3D) => {
        if (!this.supported || !this.props.enabled || !this.ready) return false;

        const { vp, r } = this.transform(s);
        const { near, far, projection } = this;

        const z = vp[2] + r;
        if (-z < near) return false;

        const { aabb, w, h, lod } = this.project(vp, r);
        if (lod >= this.levelData.length || lod < this.props.minLevel) return false;

        const { offset, size } = this.levelData[lod];

        const u = aabb[0] + w / 2;
        const v = aabb[1] + h / 2;

        const ts = size[0];
        const x = u * ts;
        const y = v * ts;

        const dx = Math.floor(x);
        const dy = Math.ceil(y);
        const dw = this.tex.getWidth();

        if (dx + 1 >= ts || dy + 1 >= ts) return false;

        const di = (ts - dy - 1) * dw + dx + offset;
        if (z > depthToViewZ(this.buffer[di], near, far, projection)) return false;

        // const di1 = (ts - dy - 1) * dw + dx + 1 + offset;
        if (z > depthToViewZ(this.buffer[di + 1], near, far, projection)) return false;

        const di2 = (ts - dy + 1 - 1) * dw + dx + offset;
        if (z > depthToViewZ(this.buffer[di2], near, far, projection)) return false;

        // const di3 = (ts - dy + 1 - 1) * dw + dx + 1 + offset;
        if (z > depthToViewZ(this.buffer[di2 + 1], near, far, projection)) return false;

        return true;
    };

    setViewport(x: number, y: number, width: number, height: number) {
        if (!this.supported) return;

        // Avoid setting dimensions to 0x0 because it causes "empty textures are not allowed" error.
        width = Math.max(width, 2);
        height = Math.max(height, 2);

        Viewport.set(this.viewport, x, y, width, height);
        const levels = Math.ceil(Math.log(Math.max(width, height)) / Math.log(2));
        if (levels === this.levelData.length) return;

        const { minLevel } = this.props;
        this.buffer = new Float32Array(Math.pow(2, levels - minLevel) * Math.pow(2, levels - 1 - minLevel));
        this.tex.define(Math.pow(2, levels - minLevel), Math.pow(2, levels - 1 - minLevel));

        for (const td of this.levelData) {
            td.framebuffer.destroy();
            td.texture.destroy();
        }

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
            if (i >= minLevel) {
                this.levelData[i].offset = offset;
                offset += td.size[0];
            }
        }

        this.clear();
    }

    setProps(props: Partial<HiZProps>) {
        if (!this.supported) return;

        if (this.props.minLevel !== props.minLevel) {
            Object.assign(this.props, props);
            const { x, y, width, height } = this.viewport;
            this.setViewport(x, y, width, height);
        } else {
            Object.assign(this.props, props);
            if (!this.props.enabled) this.clear();
        }
    }

    //

    private debug?: {
        container: HTMLDivElement
        canvas: HTMLCanvasElement
        ctx: CanvasRenderingContext2D
        rect: HTMLDivElement
    };

    private initDebug(element: HTMLElement) {
        if (!element.parentElement) return;

        const container = document.createElement('div');
        Object.assign(container.style, {
            display: 'block',
            position: 'absolute',
            pointerEvents: 'none',
        });
        element.parentElement.appendChild(container);

        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        if (!ctx) throw new Error('Could not create canvas 2d context');
        Object.assign(canvas.style, {
            width: '100%',
            height: '100%',
            imageRendering: 'pixelated',
            position: 'relative',
            pointerEvents: 'none',
        });
        container.appendChild(canvas);

        const rect = document.createElement('div');
        Object.assign(rect.style, {
            display: 'none',
            position: 'absolute',
            pointerEvents: 'none',
        });
        element.parentElement.appendChild(rect);

        this.debug = { container, canvas, ctx, rect };
    }

    private canDebug(debug: HiZPass['debug']): debug is NonNullable<HiZPass['debug']> {
        return this.supported && this.props.enabled && this.ready && !!this.debug;
    }

    private showRect(p: Vec4, occluded: boolean) {
        if (!this.canDebug(this.debug)) return;

        const { gl: { drawingBufferHeight }, pixelRatio } = this.webgl;
        const { viewport: { x, y, width, height } } = this;
        const minx = (p[0] * width + x) / pixelRatio;
        const miny = (p[1] * height - y) / pixelRatio;
        const maxx = (p[2] * width + x) / pixelRatio;
        const maxy = (p[3] * height - y) / pixelRatio;

        const oy = (drawingBufferHeight - height) / pixelRatio;

        Object.assign(this.debug.rect.style, {
            border: occluded ? 'solid red' : 'solid green',
            display: 'block',
            left: `${minx}px`,
            top: `${miny + oy}px`,
            width: `${maxx - minx}px`,
            height: `${maxy - miny}px`,
        });
    }

    private showBuffer(lod: number) {
        if (!this.canDebug(this.debug)) return;
        if (lod >= this.levelData.length || lod < this.props.minLevel) {
            this.debug.container.style.display = 'none';
            return;
        }

        const { offset, size: [tw, th] } = this.levelData[lod];
        const dw = this.tex.getWidth();

        const data = new Uint8ClampedArray(tw * th * 4);
        data.fill(255);
        for (let y = 0; y < th; ++y) {
            for (let x = 0; x < tw; ++x) {
                const i = (th - y - 1) * tw + x;
                const v = this.buffer[y * dw + x + offset] * 255;
                data[i * 4 + 0] = v;
                data[i * 4 + 3] = 255 - v;
            }
        }
        const imageData = new ImageData(data, tw, th);

        this.debug.canvas.width = imageData.width;
        this.debug.canvas.height = imageData.height;
        this.debug.ctx.putImageData(imageData, 0, 0);

        const { viewport: { x, y, width, height }, webgl: { pixelRatio } } = this;
        Object.assign(this.debug.container.style, {
            display: 'block',
            bottom: `${y / pixelRatio}px`,
            left: `${x / pixelRatio}px`,
            width: `${width / pixelRatio}px`,
            height: `${height / pixelRatio}px`,
        });
    }

    debugOcclusion(s: Sphere3D | undefined) {
        if (!this.canDebug(this.debug)) return;

        if (!s) {
            this.debug.rect.style.display = 'none';
            this.debug.container.style.display = 'none';
            return;
        }

        const occluded = this.isOccluded(s);
        const { vp, r } = this.transform(s);
        const { aabb, lod } = this.project(vp, r);

        this.showRect(aabb, occluded);
        this.showBuffer(lod);
    }

    //

    dispose() {
        if (!this.supported) return;

        this.clear();

        this.fb.destroy();
        this.tex.destroy();
        this.webgl.gl.deleteBuffer(this.buf);
        this.renderable.dispose();

        for (const td of this.levelData) {
            td.framebuffer.destroy();
            td.texture.destroy();
        }
    }

    constructor(private webgl: WebGLContext, private drawPass: DrawPass, canvas: HTMLCanvasElement | undefined, props: Partial<HiZProps>) {
        const { gl, extensions } = webgl;
        if (!isWebGL2(gl) || !extensions.colorBufferFloat) {
            if (isDebugMode) {
                console.log('Missing webgl2 and/or colorBufferFloat support required for "Hi-Z"');
            }
            this.supported = false;
            return;
        }

        this.fb = webgl.resources.framebuffer();
        this.tex = webgl.resources.texture('image-float32', 'alpha', 'float', 'nearest');

        // check red/float reading support
        this.tex.attachFramebuffer(this.fb, 0);
        const implFormat = gl.getParameter(gl.IMPLEMENTATION_COLOR_READ_FORMAT);
        const implType = gl.getParameter(gl.IMPLEMENTATION_COLOR_READ_TYPE);
        if (implFormat !== gl.RED || implType !== gl.FLOAT) {
            if (isDebugMode) {
                console.log('Missing red/float reading support required for "Hi-Z"');
            }
            this.supported = false;
            return;
        }

        this.supported = true;
        this.props = { ...PD.getDefaultValues(HiZParams), ...props };
        this.buf = getBuffer(gl);
        this.renderable = createHiZRenderable(webgl, this.drawPass.depthTextureOpaque);

        if (isDebugMode && canvas) {
            this.initDebug(canvas);
        }
    }
}
