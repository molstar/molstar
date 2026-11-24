/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderBuilder, type PassBuilder } from './pass-builder.js';
import {
  type HiZRenderable,
  createHiZRenderable,
  type HiZOccluderRenderable,
  createHiZOccluderRenderable,
} from '@molstar/graphics/gl/compute/hi-z';
import { type OccluderGrid, usesOccluderGrid } from '@molstar/graphics/gl/occluder-grid';
import { Clip } from '@molstar/core/util/clip';
import {
  ClipStorageVec4s,
  type ClipStorageValues,
  getStorageClipVersion,
  packStorageClipObjects,
} from '@molstar/graphics/gl/clip-storage';
import { projectSphere } from '@molstar/graphics/gl/project-sphere';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { fasterLog2 as _fasterLog2 } from '@molstar/core/math/approx';
import type { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec4 } from '@molstar/core/math/linear-algebra';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { Camera, ICamera } from '../camera.js';
import { StereoCamera } from '../camera/stereo.js';
import { Viewport } from '../camera/util.js';
import type { DrawPass } from './draw.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { TextureReader } from '@molstar/graphics/gl/texture-reader';
import type { OcclusionPyramid } from '@molstar/graphics/gl/compute/cell-cull';
import { createRenderPassAttachments, type RenderPassAttachments } from '@molstar/graphics/gl/render-pass-attachments';

/**
 * Where a stereo camera's right eye begins, in [0, 1] of `viewport`'s width: both eyes enabled, side by side, each of
 * the whole height, together the viewport; undefined otherwise (XR's views laid out another way)
 */
export function getStereoSplit(camera: StereoCamera, viewport: Viewport) {
  const { left, right } = camera;
  if (left.disabled || right.disabled) return undefined;
  const l = left.viewport;
  const r = right.viewport;
  if (l.x !== viewport.x || r.x !== l.x + l.width || r.x + r.width !== viewport.x + viewport.width) return undefined;
  if (l.y !== viewport.y || r.y !== viewport.y || l.height !== viewport.height || r.height !== viewport.height) {
    return undefined;
  }
  return l.width / viewport.width;
}

/** Whether `pyramid` is of `camera`'s eyes as they are (`OcclusionPyramid.second` the right eye's) */
export function isStereoPyramidOf(pyramid: OcclusionPyramid, camera: StereoCamera) {
  const { second } = pyramid;
  if (!second) return false;
  const { left, right } = camera;
  return (
    Mat4.areEqual(pyramid.view, left.view, 0) &&
    Mat4.areEqual(pyramid.projection, left.projection, 0) &&
    Mat4.areEqual(second.view, right.view, 0) &&
    Mat4.areEqual(second.projection, right.projection, 0)
  );
}

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const fasterLog2 = _fasterLog2;

// depth is reversed: the near plane at depth 1, the far plane at 0
function perspectiveDepthToViewZ(depth: number, near: number, far: number) {
  return -(near * far) / (near + (far - near) * depth);
}

function orthographicDepthToViewZ(depth: number, near: number, far: number) {
  return depth * (far - near) - far;
}

/** `rows` (`[start, end)` pairs, in any order) sorted, those that overlap or touch merged */
function mergeRows(rows: readonly number[]) {
  const pairs: [number, number][] = [];
  for (let k = 0; k < rows.length; k += 2) pairs.push([rows[k], rows[k + 1]]);
  pairs.sort((a, b) => a[0] - b[0]);
  const merged: number[] = [];
  for (const [start, end] of pairs) {
    const n = merged.length;
    if (n > 0 && start <= merged[n - 1]) merged[n - 1] = Math.max(merged[n - 1], end);
    else merged.push(start, end);
  }
  return merged;
}

function depthToViewZ(depth: number, near: number, far: number, projection: Mat4) {
  return projection[11] === -1 ? perspectiveDepthToViewZ(depth, near, far) : orthographicDepthToViewZ(depth, near, far);
}

type LevelData = {
  texture: Texture;
  target: RenderPassAttachments;
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
  twoPhase: PD.Boolean(true, {
    description:
      "GPU occlusion culling catching up within the frame: what the previous frame's depth left out is tested against the frame's depth so far and drawn where visible. Exact while the camera moves; off, what comes into view shows a frame late, or not until the next draw.",
  }),
  occluderOpticalDepth: PD.Numeric(
    5,
    { min: 1, max: 20, step: 1 },
    {
      label: 'Occluder Optical Depth',
      description:
        "Where content with an occluder grid (particles' densities) hides what is behind it, by the light passing: 37% at 1, 5% at 3, 0.1% at 7. Lower culls more behind thin content (faster, sometimes a gap), higher only behind dense. WebGPU.",
    },
  ),
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
  /** the second view of the levels being built (`OcclusionPyramid.second`): its split 1 without one */
  private buildSplit = 1;
  private readonly buildSecondView = Mat4();
  private readonly buildSecondProjection = Mat4();
  private pyramid?: OcclusionPyramid;
  /** the builds of `pyramid` so far (`OcclusionPyramid.version`) */
  private pyramidVersion = 0;

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
  /** the camera of the last frame culled by the CPU test (`noteCulled`): its view, projection scale, near and far */
  private readonly culledView = Mat4();
  private readonly culledLens = Vec4();
  /** that frame culled against levels of another view: what came into view hidden by them (`noteCulled`) */
  private culledStale = false;
  /** a frame of that camera drawn again once (`takeCatchUp`): not again until the camera moves */
  private caughtUp = false;
  private catchUp = false;

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

  /**
   * `complete`: no draw of the frame's depth skipped (e.g. a pipeline compiling): built also while other work is pending
   * (a readback: a pick, the CPU test's own), else only once nothing is; `read`: for the CPU test as well (`enabled`).
   * A stereo camera's eyes side by side over the viewport: one pyramid of both (`OcclusionPyramid.second`), for the
   * tests on the GPU only (not read); other layouts of its eyes build none
   */
  render(camera: Camera | StereoCamera, complete = false, read = true) {
    // unsupported (e.g. WebGL1: no float targets or asynchronous reads): no reader, nothing to build
    if (!this.supported) return;
    if (!read) {
      // the CPU test's levels (read back) of an earlier frame no more: frames of levels not read (stereo's eyes,
      // illumination) test nothing on the CPU, and the next frame that does (stereo or illumination off) none until
      // its own are read
      this.reader.reset();
      this.reading = false;
      this.frameLag = 0;
      this.ready = false;
    }
    this.plan(this.renderBuilder, camera, this.props, this.gpuCulling, complete, read);
  }

  /**
   * The levels of the opaque depth drawn so far with `camera`, planned within the draw pass (`b`), for the second phase
   * of GPU culling (`RenderContext.occlusion`); undefined where not built (a draw of the depth skipped, `gpu` off). The
   * pass builds them again after the frame, of its whole depth: the next frame's and the CPU test's
   */
  planFrame(b: PassBuilder, camera: ICamera, complete: boolean): OcclusionPyramid | undefined {
    if (!complete || !this.props.gpu) return undefined;
    const last = this.pyramid;
    this.plan(b, camera, this.props, true, true, false);
    return this.pyramid === last ? undefined : this.pyramid;
  }

  /** the occluder grid the first level is primed from, its texture (WebGPU; made at the first grid) */
  private occluder?: OccluderGrid;
  private occluderTexture?: Texture;
  private occluderUploaded?: OccluderGrid;
  /** the rows of the atlas an update wrote (`[start, end)` pairs), a scratch */
  private occluderRows: number[] = [];
  /** the slices' atlas of the last grid uploaded: of the next where the same size */
  private occluderAtlas?: Float32Array;

  private occluderRenderable?: HiZOccluderRenderable;
  private readonly invViewProjection = Mat4();
  private readonly viewProjection = Mat4();

  /**
   * What the frame's content hides (`OccluderGrid`, a visible object's): the first level primed from it, a texel's
   * depth the nearer of the drawn depth's and where its ray through the grid reaches the occlusion's threshold
   * (`HiZProps.occluderOpticalDepth`), so the levels hold through the gaps between molecules. WebGPU; none: as drawn
   */
  setOccluder(grid: OccluderGrid | undefined, clip?: ClipStorageValues) {
    // WebGPU's program primes the level (WebGL's builds it plain: no grid uploaded there)
    this.occluder = this.supported && usesOccluderGrid(this.graphics) ? grid : undefined;
    // the grid's object's clip objects (its values: a drag's moved cut too, its uniforms the same frame)
    this.occluderClip = this.occluder ? clip : undefined;
  }

  /** the clip values of the occluder's object (`setOccluder`), and the version of them packed last */
  private occluderClip?: ClipStorageValues;
  private occluderClipVersion = '';
  private occluderClipPackedOf?: ClipStorageValues;
  private readonly occluderClipPacked = new Float32Array(Clip.MaxObjectCount * ClipStorageVec4s * 4);

  /** the occluder renderable for the first level with `grid` uploaded, made where none */
  private getOccluderRenderable(grid: OccluderGrid) {
    const [nx, ny, nz] = grid.dimensions;
    const columns = Math.ceil(Math.sqrt(nz));
    if (!this.occluderTexture)
      this.occluderTexture = this.graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
    if (this.occluderUploaded !== grid) {
      // the slices along z side by side, `columns` a row
      const rows = Math.ceil(nz / columns);
      const width = columns * nx;
      const height = rows * ny;
      // an update of the grid uploaded last (its array written in place): its cells' rows alone (a view's step: a few
      // tiles' of millions of cells), not the atlas again
      const last = this.occluderUploaded;
      const atlas0 = this.occluderAtlas;
      if (
        grid.update &&
        last?.density === grid.density &&
        last.version === grid.update.from &&
        atlas0?.length === width * height &&
        this.occluderTexture.loadRows
      ) {
        const boxes = grid.update.boxes;
        const written = this.occluderRows;
        written.length = 0;
        for (let k = 0; k < boxes.length; k += 6) {
          const [x0, y0, z0, x1, y1, z1] = [
            boxes[k],
            boxes[k + 1],
            boxes[k + 2],
            boxes[k + 3],
            boxes[k + 4],
            boxes[k + 5],
          ];
          for (let z = z0; z < z1; ++z) {
            const ox = (z % columns) * nx;
            const oy = Math.floor(z / columns) * ny;
            for (let y = y0; y < y1; ++y) {
              const from = (z * ny + y) * nx;
              atlas0.set(grid.density.subarray(from + x0, from + x1), (oy + y) * width + ox + x0);
            }
            written.push(oy + y0, oy + y1);
          }
        }
        this.occluderTexture.loadRows({ array: atlas0, width, height }, mergeRows(written));
        this.occluderUploaded = grid;
      }
    }
    if (this.occluderUploaded !== grid) {
      const rows = Math.ceil(nz / columns);
      const width = columns * nx;
      const height = rows * ny;
      // the last's array where of the same size (a view's steps: a grid each, of millions of cells)
      const atlas =
        this.occluderAtlas?.length === width * height
          ? this.occluderAtlas.fill(0)
          : (this.occluderAtlas = new Float32Array(width * height));

      for (let z = 0; z < nz; ++z) {
        const ox = (z % columns) * nx;
        const oy = Math.floor(z / columns) * ny;
        for (let y = 0; y < ny; ++y) {
          atlas.set(grid.density.subarray((z * ny + y) * nx, (z * ny + y + 1) * nx), (oy + y) * width + ox);
        }
      }
      this.occluderTexture.load({ array: atlas, width, height });
      this.occluderUploaded = grid;
    }
    if (!this.occluderRenderable) {
      this.occluderRenderable = createHiZOccluderRenderable(
        this.graphics,
        this.drawPass.depthTextureOpaque,
        this.occluderTexture,
      );
    }
    const v = this.occluderRenderable.values;
    ValueCell.update(v.uGridMin, Vec3.copy(v.uGridMin.ref.value, grid.min));
    ValueCell.update(v.uGridMax, Vec3.copy(v.uGridMax.ref.value, grid.max));
    ValueCell.update(v.uGridDims, Vec3.set(v.uGridDims.ref.value, nx, ny, nz));
    ValueCell.updateIfChanged(v.uGridColumns, columns);
    ValueCell.updateIfChanged(v.uOpticalDepth, this.props.occluderOpticalDepth);
    // its clip objects packed where they changed (a cut moved: a frame's few floats)
    const clip = this.occluderClip;
    // its versions of the object packed last (another object's cells can be of the same versions)
    const version = clip ? getStorageClipVersion(clip) : '';
    if (clip !== this.occluderClipPackedOf || version !== this.occluderClipVersion) {
      this.occluderClipPackedOf = clip;
      this.occluderClipVersion = version;
      const data = v.uClipData.ref.value;
      const packed = this.occluderClipPacked;
      const count = clip ? packStorageClipObjects(clip, packed) : 0;
      for (let k = 0; k < data.length; ++k) data[k] = packed[k];
      ValueCell.update(v.uClipData, data);
      ValueCell.updateIfChanged(v.uClipCount, count);
    }
    return this.occluderRenderable;
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
  plan(
    b: PassBuilder,
    camera: ICamera | StereoCamera,
    props: HiZProps = this.props,
    gpuCulling = this.gpuCulling,
    complete = false,
    /** the CPU test's read of these levels: not of a frame's depth so far (`planFrame`) */
    read = true,
  ) {
    if (!this.isActive(props, gpuCulling)) return;
    // no viewport yet (`setViewport`): no levels to build or read
    if (this.levelData.length <= this.props.minLevel) return;
    const { pending } = this.graphics.resources;
    if (b.drawing) {
      // draws are skipped while work is pending (e.g. a pipeline compiling), levels missing them would occlude
      // everything; a depth known complete is built also while other work is (a readback)
      if (!complete && !pending.isReady()) return;
      const stereo = camera instanceof StereoCamera;
      const split = stereo ? getStereoSplit(camera, this.viewport) : 1;
      if (split === undefined) {
        this.clear();
        return;
      }
      // built every frame for the tests on the GPU (`gpuPyramid`), also while the last read for the CPU test is pending
      const first = stereo ? camera.left : camera;
      this.buildNear = first.near;
      this.buildFar = first.far;
      Mat4.copy(this.buildView, first.view);
      Mat4.copy(this.buildProjection, first.projection);
      this.buildSplit = split;
      if (stereo) {
        Mat4.copy(this.buildSecondView, camera.right.view);
        Mat4.copy(this.buildSecondProjection, camera.right.projection);
        // the CPU test reads one view's levels
        read = false;
      }
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
    // the first level primed from the occluder grid where there is one (one view: not a stereo pair's)
    const occluder = this.occluder && this.buildSplit >= 1 ? this.getOccluderRenderable(this.occluder) : undefined;
    if (occluder) {
      const o = occluder.values;
      Mat4.mul(this.viewProjection, this.buildProjection, this.buildView);
      Mat4.invert(this.invViewProjection, this.viewProjection);
      ValueCell.update(o.uViewProjection, Mat4.copy(o.uViewProjection.ref.value, this.viewProjection));
      ValueCell.update(o.uInvViewProjection, Mat4.copy(o.uInvViewProjection.ref.value, this.invViewProjection));
      ValueCell.update(
        o.uViewportRect,
        Vec4.set(
          o.uViewportRect.ref.value,
          this.viewport.x / drawingBuffer.width,
          this.viewport.y / drawingBuffer.height,
          this.viewport.width / drawingBuffer.width,
          this.viewport.height / drawingBuffer.height,
        ),
      );
    }
    for (let i = 0; i < levels; ++i) {
      const td = this.levelData[i];
      b.begin(td.target, { label: 'hi-z' });
      // the first level the occluder's where primed: its values as the level's
      const lv = i === 0 && occluder ? occluder.values : v;
      if (lv !== v) {
        ValueCell.update(lv.uTexSize, td.size);
        ValueCell.update(lv.uInvSize, Vec2.set(this.levelZeroInvSize, 1 / s, 1 / s));
        ValueCell.update(
          lv.uRect,
          Vec4.set(
            lv.uRect.ref.value,
            this.viewport.x / drawingBuffer.width,
            this.viewport.y / drawingBuffer.height,
            td.size[0] / s,
            td.size[1] / s,
          ),
        );
        ValueCell.update(lv.tPreviousLevel, this.drawPass.depthTextureOpaque);
        if (b.drawing) {
          graphics.setViewport(0, 0, td.size[0], td.size[1]);
          graphics.clear({ color: [0, 0, 0, 0] });
        }
        occluder!.update();
        b.quad(occluder!, RenderState.Default);
        if (b.drawing && i >= this.props.minLevel) {
          this.tex.copyFromTexture(td.texture, 0, 0, td.size[0], td.size[1], td.offset, 0);
        }
        continue;
      }
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
        this.tex.copyFromTexture(td.texture, 0, 0, td.size[0], td.size[1], td.offset, 0);
      }
    }

    //

    if (b.drawing && pending.trackedCount === trackedCount) {
      const { width, height } = this.viewport;
      const last = this.pyramid?.second;
      const second =
        this.buildSplit < 1
          ? {
              view: Mat4.copy(last?.view ?? Mat4(), this.buildSecondView),
              projection: Mat4.copy(last?.projection ?? Mat4(), this.buildSecondProjection),
              split: this.buildSplit,
            }
          : undefined;
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
        second,
        version: ++this.pyramidVersion,
      };
      if (read && this.props.enabled && !this.reading) {
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
      // the last frame culled against older levels (a camera stopped after moving): drawn again with these, once
      if (this.culledStale) {
        this.culledStale = false;
        this.caughtUp = true;
        this.catchUp = true;
      }
    }
  }

  /**
   * After a frame culled by the CPU test (`occlusionTest`) with `camera`: whether its levels were of another view. The
   * test lags behind (levels read back a frame or more later, no second phase), so what came into view is missing in
   * that frame; with the camera stopped no frame draws it, so the next levels to arrive draw one (`takeCatchUp`), once
   * a camera pose. The view compared by its matrix and the projection's scale (not a multi-sample's jitter)
   */
  noteCulled(camera: Camera) {
    if (!this.supported || !this.props.enabled) return;
    const { view, projection, near, far } = camera;
    const moved =
      !Mat4.areEqual(this.culledView, view, 0) ||
      this.culledLens[0] !== projection[0] ||
      this.culledLens[1] !== projection[5] ||
      this.culledLens[2] !== near ||
      this.culledLens[3] !== far;
    if (moved) {
      Mat4.copy(this.culledView, view);
      Vec4.set(this.culledLens, projection[0], projection[5], near, far);
      this.caughtUp = false;
    }
    const stale =
      this.ready &&
      (!Mat4.areEqual(this.view, view, 0) ||
        this.projection[0] !== projection[0] ||
        this.projection[5] !== projection[5] ||
        this.near !== near ||
        this.far !== far);
    this.culledStale = stale && !this.caughtUp;
  }

  /** whether a frame is to be drawn again with the levels that arrived (`noteCulled`); taken once */
  takeCatchUp() {
    const value = this.catchUp;
    this.catchUp = false;
    return value;
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

    const projected = this.project(vp, r);
    const { aabb } = projected;
    // below the lowest level kept, tested there (its texels larger: as conservative)
    let lod = Math.max(projected.lod, this.props.minLevel);

    // the texels the rectangle covers (rows top to bottom, as the rectangle's y): at most two per axis, a level up
    // where it straddles three (its size up to two texels at the level of its size); the 2x2 texels around its
    // center miss a side where it is off the texel's center
    let ts = 0;
    let x0 = 0;
    let x1 = 0;
    let y0 = 0;
    let y1 = 0;
    for (; ; ++lod) {
      if (lod >= this.levelData.length) return false;
      ts = this.levelData[lod].size[0];
      x0 = Math.floor(aabb[0] * ts);
      x1 = Math.floor(aabb[2] * ts);
      y0 = Math.floor(aabb[1] * ts);
      y1 = Math.floor(aabb[3] * ts);
      if (x1 - x0 <= 1 && y1 - y0 <= 1) break;
    }
    if (x0 < 0 || y0 < 0 || x1 >= ts || y1 >= ts) return false;
    const { offset } = this.levelData[lod];

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
      td.target.destroy();
      td.texture.destroy();
    }

    this.levelData.length = 0;
    for (let i = 0; i < levels; ++i) {
      const levelSize = Math.pow(2, levels - i - 1);
      const size = Vec2.create(levelSize, levelSize);
      const invSize = Vec2.create(1 / levelSize, 1 / levelSize);
      const texture = this.graphics.resources.texture('image-float32', 'r', 'float', 'nearest');
      texture.define(levelSize, levelSize);
      const target = createRenderPassAttachments(this.graphics.resources, [texture]);
      this.levelData.push({ texture, target, size, invSize, offset: 0 });
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
      td.target.destroy();
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
