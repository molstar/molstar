/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderBuilder, type PassBuilder } from './pass-builder.js';
import { type HiZRenderable, createHiZRenderable } from '@molstar/graphics/gl/compute/hi-z';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import type { Texture } from '@molstar/graphics/gl/texture';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { fasterLog2 as _fasterLog2 } from '@molstar/core/math/approx';
import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec4 } from '@molstar/core/math/linear-algebra';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { Camera } from '../camera.js';
import { Viewport } from '../camera/util.js';
import type { DrawPass } from './draw.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { TextureReader } from '@molstar/graphics/gl/texture-reader';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const v4set = Vec4.set;
const fasterLog2 = _fasterLog2;

// depth is reversed: the near plane at depth 1, the far plane at 0
function perspectiveDepthToViewZ(depth: number, near: number, far: number) {
  return -(near * far) / (near + (far - near) * depth);
}

function orthographicDepthToViewZ(depth: number, near: number, far: number) {
  return depth * (far - near) - far;
}

function depthToViewZ(depth: number, near: number, far: number, projection: Mat4) {
  return projection[11] === -1 ? perspectiveDepthToViewZ(depth, near, far) : orthographicDepthToViewZ(depth, near, far);
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

  return v4set(out, maxx * -0.5 + 0.5, miny * 0.5 + 0.5, minx * -0.5 + 0.5, maxy * 0.5 + 0.5);
}

function orthographicProjectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
  const sx = projection[0];
  const sy = projection[5];

  const minx = (p[0] + r) * sx;
  const maxx = (p[0] - r) * sx;
  const miny = (p[1] + r) * sy;
  const maxy = (p[1] - r) * sy;

  return v4set(out, maxx * 0.5 + 0.5, miny * -0.5 + 0.5, minx * 0.5 + 0.5, maxy * -0.5 + 0.5);
}

function projectSphere(out: Vec4, p: Vec3, r: number, projection: Mat4) {
  return projection[11] === -1
    ? perspectiveProjectSphere(out, p, r, projection)
    : orthographicProjectSphere(out, p, r, projection);
}

type LevelData = {
  texture: Texture;
  framebuffer: Framebuffer;
  size: Vec2;
  invSize: Vec2;
  offset: number;
}[];

export const HiZParams = {
  enabled: PD.Boolean(false, {
    description: 'Hierarchical Z-buffer occlusion culling. Needs asynchronous texture reads (WebGL2 or WebGPU).',
  }),
  maxFrameLag: PD.Numeric(
    10,
    { min: 1, max: 30, step: 1 },
    { description: 'Maximum number of frames to wait for Z-buffer data.' },
  ),
  minLevel: PD.Numeric(3, { min: 1, max: 10, step: 1 }),
  gpu: PD.Boolean(true, {
    description:
      "Occlusion culling on the GPU of the objects culled there (gpuCulling), against the previous frame's depth without reading it back; also with Hi-Z off.",
  }),
};
export type HiZProps = PD.Values<typeof HiZParams>;

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

  /** the camera of the levels being built and of those last built (`gpuPyramid`) */
  private buildNear = 0;
  private buildFar = 0;
  private readonly buildView = Mat4();
  private readonly buildProjection = Mat4();
  private pyramid?: OcclusionPyramid;

  private readonly aabb = Vec4();
  private readonly levelZeroInvSize = Vec2();
  private readonly vp = Vec3();

  private readonly levelData: LevelData = [];
  /** `minLevel` the `levelData` offsets and `tex` size were computed for */
  private layoutMinLevel = -1;
  private readonly reader: TextureReader<Float32Array>;
  private readonly tex: Texture;
  private readonly renderable: HiZRenderable;
  private readonly renderBuilder: PassBuilder;
  private readonly supported: boolean;

  /** whether a read of `tex` was started and is pending */
  private reading = false;
  /** whether objects are culled on the GPU (`setGpuCulling`): the levels built for them also with Hi-Z off */
  private gpuCulling = false;
  private buffer = new Float32Array(0);
  private frameLag = 0;
  private ready = false;

  readonly props: HiZProps;

  getByteCount() {
    if (!this.supported) return 0;
    return (
      this.tex.getByteCount() +
      this.reader.getByteCount() +
      this.levelData.reduce((sum, l) => sum + l.texture.getByteCount(), 0)
    );
  }

  clear() {
    if (!this.supported) return;

    this.reader.reset();
    this.reading = false;
    this.frameLag = 0;
    this.ready = false;
    this.pyramid = undefined;

    if (this.debug) {
      this.debug.rect.style.display = 'none';
      this.debug.container.style.display = 'none';
    }
  }

  render(camera: Camera) {
    this.plan(this.renderBuilder, camera);
  }

  /** Whether objects are culled on the GPU this frame (`Renderer.gpuCulledCount`), see `gpu` */
  setGpuCulling(value: boolean) {
    this.gpuCulling = value;
  }

  /** Built: for the CPU test (`enabled`), or for the GPU's alone (`gpu` with objects culled there) */
  private get active() {
    return this.isActive(this.props, this.gpuCulling);
  }

  private isActive(props: HiZProps, gpuCulling: boolean) {
    return this.supported && (props.enabled || (props.gpu && gpuCulling));
  }

  /**
   * What `render` draws, see `PassBuilder`; preparing covers the first level (from the scene depth) and one of the
   * following ones (each from the level before), and keeps the reading state. `props` and `gpuCulling` other than
   * the pass's for preparing what a canvas is not drawing yet (`Canvas3D.prepare`, `precompile`: objects culled on
   * the GPU before any is drawn)
   */
  plan(b: PassBuilder, camera: Camera, props: HiZProps = this.props, gpuCulling = this.gpuCulling) {
    if (!this.isActive(props, gpuCulling)) return;
    // no viewport yet (`setViewport`): no levels to build or read
    if (this.levelData.length <= this.props.minLevel) return;
    const { pending } = this.graphics.resources;
    if (b.drawing) {
      // draws are skipped while work is pending (e.g. a pipeline compiling), levels missing them would occlude
      // everything
      if (!pending.isReady()) return;
      // built every frame for the tests on the GPU (`gpuPyramid`), also while the last read for the CPU test is pending
      this.buildNear = camera.near;
      this.buildFar = camera.far;
      Mat4.copy(this.buildView, camera.view);
      Mat4.copy(this.buildProjection, camera.projection);
    }

    const { graphics } = this;
    const timing = isTimingMode && b.drawing;
    if (timing) graphics.timer.group('hi-Z');

    if (b.drawing) {
      graphics.setRenderState(RenderState.Default);
      graphics.disableScissor();
    }
    const trackedCount = pending.trackedCount;

    //

    const v = this.renderable.values;
    const drawingBuffer = graphics.getDrawingBufferSize();
    const s = Math.pow(2, Math.ceil(Math.log(Math.max(drawingBuffer.width, drawingBuffer.height)) / Math.log(2)) - 1);

    const levels = b.drawing ? this.levelData.length : Math.min(this.levelData.length, 2);
    for (let i = 0; i < levels; ++i) {
      const td = this.levelData[i];
      b.begin(td.framebuffer, { label: 'hi-z' });
      ValueCell.update(v.uTexSize, td.size);

      if (i > 0) {
        ValueCell.update(v.uInvSize, td.invSize);
        ValueCell.update(v.uRect, Vec4.set(v.uRect.ref.value, 0, 0, 1, 1));
        ValueCell.update(v.tPreviousLevel, this.levelData[i - 1].texture);
      } else {
        // a vector of its own: the cell holds a level's `invSize` after that level, set in place it would change it
        ValueCell.update(v.uInvSize, Vec2.set(this.levelZeroInvSize, 1 / s, 1 / s));
        // the level covers `td.size / s` of the drawing buffer from the top-left corner of the viewport
        ValueCell.update(
          v.uRect,
          Vec4.set(
            v.uRect.ref.value,
            this.viewport.x / drawingBuffer.width,
            this.viewport.y / drawingBuffer.height,
            td.size[0] / s,
            td.size[1] / s,
          ),
        );
        ValueCell.update(v.tPreviousLevel, this.drawPass.depthTextureOpaque);
      }

      if (b.drawing) {
        graphics.setViewport(0, 0, td.size[0], td.size[1]);
        graphics.clear({ color: [0, 0, 0, 0] });
      }
      this.renderable.update();
      b.quad(this.renderable, RenderState.Default);

      if (b.drawing && i >= this.props.minLevel) {
        this.tex.copyFromFramebuffer(td.framebuffer, 0, 0, td.size[0], td.size[1], td.offset, 0);
      }
    }

    //

    if (b.drawing && pending.trackedCount === trackedCount) {
      const { width, height } = this.viewport;
      this.pyramid = {
        texture: this.tex,
        levelCount: this.levelData.length,
        minLevel: this.props.minLevel,
        width,
        height,
        view: Mat4.copy(this.pyramid?.view ?? Mat4(), this.buildView),
        projection: Mat4.copy(this.pyramid?.projection ?? Mat4(), this.buildProjection),
        near: this.buildNear,
        far: this.buildFar,
      };
      if (this.props.enabled && !this.reading) {
        // the CPU test reads these levels with this camera once read back (`tick`)
        this.nextNear = this.buildNear;
        this.nextFar = this.buildFar;
        Mat4.copy(this.nextView, this.buildView);
        Mat4.copy(this.nextProjection, this.buildProjection);
        this.reader.start([this.tex], 0, 0, this.tex.getWidth(), this.tex.getHeight());
        this.reading = true;
        this.frameLag = 0;
      }
    }

    if (timing) graphics.timer.groupEnd('hi-Z');
  }

  /**
   * The occlusion test for the renderer (`Renderer.setOcclusionTest`), null while Hi-Z is off or unsupported: no test
   * then, so that what culling read stays unchanged (e.g. draw lists replay without planning)
   */
  get occlusionTest() {
    return this.supported && this.props.enabled ? this.isOccluded : null;
  }

  /** The levels last built, for tests on the GPU (`OcclusionPyramid`): undefined while not built (see `active`) */
  get gpuPyramid(): OcclusionPyramid | undefined {
    return this.active ? this.pyramid : undefined;
  }

  /** A readback under way: `tick` polls it. */
  get isPending() {
    return this.supported && this.props.enabled && this.reading;
  }

  tick() {
    if (!this.supported || !this.props.enabled || !this.reading) return;

    const status = this.reader.check([this.buffer]);
    if (status === 'failed' || (status === 'pending' && this.frameLag >= this.props.maxFrameLag)) {
      // console.log(`failed to get buffer data after ${this.frameLag + 1} frames`);
      this.reader.reset();
      this.reading = false;
      this.frameLag = 0;
      this.ready = false;
    } else if (status === 'pending') {
      this.frameLag += 1;
      // console.log(`waiting for buffer data for ${this.frameLag} frames`);
    } else {
      // console.log(`got buffer data after ${this.frameLag + 1} frames`);
      this.reading = false;
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
    const r = s.radius * 1.2 + 1.52;
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

    const { aabb, lod } = this.project(vp, r);
    if (lod >= this.levelData.length || lod < this.props.minLevel) return false;

    const { offset, size } = this.levelData[lod];
    const ts = size[0];

    // the texels the rectangle covers (rows top to bottom, as the rectangle's y): at most two per axis at this level,
    // its size at most a texel; the 2x2 texels around its center miss a side where it is off the texel's center
    const x0 = Math.floor(aabb[0] * ts);
    const x1 = Math.floor(aabb[2] * ts);
    const y0 = Math.floor(aabb[1] * ts);
    const y1 = Math.floor(aabb[3] * ts);
    if (x0 < 0 || y0 < 0 || x1 >= ts || y1 >= ts || x1 - x0 > 1 || y1 - y0 > 1) return false;

    // the levels side by side at the top of the texture
    const dw = this.tex.getWidth();
    const { buffer } = this;
    for (const y of [y0, y1]) {
      for (const x of [x0, x1]) {
        if (z > depthToViewZ(buffer[y * dw + x + offset], near, far, projection)) return false;
      }
    }
    return true;
  };

  setViewport(x: number, y: number, width: number, height: number) {
    if (!this.supported) return;

    // Avoid setting dimensions to 0x0 because it causes "empty textures are not allowed" error.
    width = Math.max(width, 2);
    height = Math.max(height, 2);

    Viewport.set(this.viewport, x, y, width, height);
    const levels = Math.ceil(Math.log(Math.max(width, height)) / Math.log(2));
    if (levels === this.levelData.length && this.layoutMinLevel === this.props.minLevel) return;

    const { minLevel } = this.props;
    this.layoutMinLevel = minLevel;
    this.buffer = new Float32Array(Math.pow(2, levels - minLevel) * Math.pow(2, levels - 1 - minLevel));
    this.tex.define(Math.pow(2, levels - minLevel), Math.pow(2, levels - 1 - minLevel));

    for (const td of this.levelData) {
      td.framebuffer.destroy();
      td.texture.destroy();
    }

    this.levelData.length = 0;
    for (let i = 0; i < levels; ++i) {
      const framebuffer = this.graphics.resources.framebuffer();
      const levelSize = Math.pow(2, levels - i - 1);
      const size = Vec2.create(levelSize, levelSize);
      const invSize = Vec2.create(1 / levelSize, 1 / levelSize);
      const texture = this.graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
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
    container: HTMLDivElement;
    canvas: HTMLCanvasElement;
    ctx: CanvasRenderingContext2D;
    rect: HTMLDivElement;
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

    const { pixelRatio } = this.graphics;
    const {
      viewport: { x, y, width, height },
    } = this;
    const minx = (p[0] * width + x) / pixelRatio;
    const miny = (p[1] * height + y) / pixelRatio;
    const maxx = (p[2] * width + x) / pixelRatio;
    const maxy = (p[3] * height + y) / pixelRatio;

    Object.assign(this.debug.rect.style, {
      border: occluded ? 'solid red' : 'solid green',
      display: 'block',
      left: `${minx}px`,
      top: `${miny}px`,
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

    const {
      offset,
      size: [tw, th],
    } = this.levelData[lod];
    const dw = this.tex.getWidth();

    const data = new Uint8ClampedArray(tw * th * 4);
    data.fill(255);
    for (let y = 0; y < th; ++y) {
      for (let x = 0; x < tw; ++x) {
        const i = y * tw + x;
        const v = this.buffer[y * dw + x + offset] * 255;
        data[i * 4 + 0] = v;
        data[i * 4 + 3] = 255 - v;
      }
    }
    const imageData = new ImageData(data, tw, th);

    this.debug.canvas.width = imageData.width;
    this.debug.canvas.height = imageData.height;
    this.debug.ctx.putImageData(imageData, 0, 0);

    const {
      viewport: { x, y, width, height },
      graphics: { pixelRatio },
    } = this;
    Object.assign(this.debug.container.style, {
      display: 'block',
      top: `${y / pixelRatio}px`,
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

    this.reader.dispose();
    this.tex.destroy();
    this.renderable.dispose();

    for (const td of this.levelData) {
      td.framebuffer.destroy();
      td.texture.destroy();
    }
  }

  constructor(
    private graphics: GraphicsContext,
    private drawPass: DrawPass,
    canvas: HTMLCanvasElement | undefined,
    props: Partial<HiZProps>,
  ) {
    this.renderBuilder = createRenderBuilder(graphics);
    this.props = { ...PD.getDefaultValues(HiZParams), ...props };

    if (!graphics.capabilities.floatRenderTarget) {
      if (isDebugMode) console.log('Missing float render targets required for "Hi-Z"');
      this.supported = false;
      return;
    }
    const reader = graphics.createTextureReader({ format: 'r32f' });
    if (!reader.asynchronous) {
      if (isDebugMode) console.log('Missing asynchronous texture reads (WebGL2) required for "Hi-Z"');
      reader.dispose();
      this.supported = false;
      return;
    }

    this.supported = true;
    this.reader = reader;
    this.tex = graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
    this.renderable = createHiZRenderable(graphics, this.drawPass.depthTextureOpaque);

    if (isDebugMode && canvas) {
      this.initDebug(canvas);
    }
  }
}
