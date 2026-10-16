/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { radToDeg } from '@molstar/core/math/misc';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Camera } from '@molstar/graphics/canvas3d/camera';
import type { Canvas3D } from '@molstar/graphics/canvas3d/canvas3d';
import { ApiArea, ApiResult } from './area.js';
import { EveryStructureParam, OptionalParam, StructureParam } from './params.js';
import { isFailure, labelOf, someStructures } from './scope.js';
import { LanguageParam, TargetParam } from './target-param.js';

// `camera` (operations.md §2, plan 2.3): view state, never a history step (sequence.md §7.3). `reset`, `frame` (a
// target), `zoom`, `orient`, `set` (position, target, up: v5's camera info, parity.md; the projection is `canvas.camera`'s) and `get`. What moves the camera needs the canvas (`modules: ['canvas']`); a headless
// viewer says so instead of running.

/** The camera in values, with its projection (the Export card's Camera section); `camera.get` gives `CameraView`. */
export interface CameraValue {
  readonly mode: 'perspective' | 'orthographic';
  /** Field of view in degrees. */
  readonly fov: number;
  /** Where the camera is, what it looks at, and its up direction (Å, rounded to 0.01). */
  readonly position: readonly [number, number, number];
  readonly target: readonly [number, number, number];
  readonly up: readonly [number, number, number];
  /** From the camera to its target. */
  readonly distance: number;
  /** The radius of what the camera shows. */
  readonly radius: number;
  /** The radius of the scene. */
  readonly sceneRadius: number;
  /** The camera shows the whole scene (best effort: the radius against the scene's). */
  readonly showsAll: boolean;
}

const ZoomDurationMs = 250;

export class CameraApi extends ApiArea {
  readonly id = 'camera';

  readonly reset = this.operation({
    id: 'camera.reset',
    title: 'Reset the camera',
    params: { animate: AnimateParam() },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ animate }, ctx) {
      ctx.plugin.managers.camera.reset(undefined, duration(ctx, animate));
      return ApiResult.ok('Camera reset.', { view: true });
    },
  });

  readonly frame = this.operation({
    id: 'camera.frame',
    title: 'Frame',
    params: {
      target: TargetParam('What to bring into view: ligand, A:120-130, resn HEM, …'),
      structure: StructureParam(),
      language: LanguageParam(),
      animate: AnimateParam(),
    },
    history: false,
    modules: ['canvas'],
    async run({ target, structure, language, animate }, ctx) {
      const resolved = await ctx.targets.resolve(target, {
        ...(structure ? { structure } : {}),
        ...(language ? { language } : {}),
      });
      if (isFailure(resolved)) return resolved;
      if (resolved.empty) return ApiResult.fail('not-found', resolved.summary, 'Nothing was changed.');
      const loci = resolved.parts.filter((p) => p.selected.atoms > 0).map((p) => p.loci);
      ctx.plugin.managers.camera.focusLoci(loci.length === 1 ? loci[0] : loci, { durationMs: duration(ctx, animate) });
      return ApiResult.ok(`Framed ${resolved.summary}.`, { view: true });
    },
  });

  readonly zoom = this.operation({
    id: 'camera.zoom',
    title: 'Zoom',
    params: {
      factor: PD.Numeric(1.5, { min: 0.1, max: 10, step: 0.1 }, { description: 'Above 1 zooms in, below 1 out' }),
      animate: AnimateParam(),
    },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ factor, animate }, ctx) {
      const canvas = ctx.plugin.canvas3d;
      if (!canvas) return ApiResult.fail('unavailable', 'There is no canvas.');
      if (factor === 1) return ApiResult.already('The zoom is unchanged.');
      // The same view direction and target, the radius it shows scaled.
      const { target, radius } = canvas.camera.state;
      const snapshot = canvas.camera.getFocus(target, Math.max(radius, 0.01) / factor);
      ctx.plugin.managers.camera.setSnapshot(snapshot, duration(ctx, animate));
      return ApiResult.ok(factor > 1 ? `Zoomed in ${factor}×.` : `Zoomed out ${round(1 / factor)}×.`, { view: true });
    },
  });

  readonly orient = this.operation({
    id: 'camera.orient',
    title: 'Orient along the principal axes',
    params: { structure: EveryStructureParam(), animate: AnimateParam() },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ structure, animate }, ctx) {
      const named = structure ? someStructures(ctx, structure) : undefined;
      if (named && isFailure(named)) return named;
      const data = named?.flatMap((s) => (s.cell.obj ? [s.cell.obj.data] : []));
      ctx.plugin.managers.camera.orientAxes(data, duration(ctx, animate));
      return ApiResult.ok(
        named?.length
          ? `Oriented the view along the axes of ${named.map((s) => labelOf(ctx, s)).join(', ')}.`
          : 'Oriented the view along the principal axes.',
        { view: true },
      );
    },
  });

  // v5's "Align Cartesian axes" (the Reset menu): the view along x, y and z again.
  readonly resetAxes = this.operation({
    id: 'camera.reset-axes',
    title: 'Align the view with the x, y and z axes',
    params: { animate: AnimateParam() },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ animate }, ctx) {
      ctx.plugin.managers.camera.resetAxes(duration(ctx, animate));
      return ApiResult.ok('Aligned the view with the x, y and z axes.', { view: true });
    },
  });

  // v5's camera info, settable: what `camera.get` gives, back (MolViewSpec's camera is position, target and up). The
  // vectors are numbers or "x, y, z" (a form's and a tool's text); what `get` only reads is taken and ignored, so
  // `set((await get()).value)` restores a camera.
  readonly set = this.operation({
    id: 'camera.set',
    title: 'Set the camera',
    params: {
      position: OptionalParam(VecParam('Where the camera is: "x, y, z" in Å')),
      target: OptionalParam(VecParam('What it looks at: "x, y, z" in Å')),
      up: OptionalParam(VecParam('Its up direction: "x, y, z"')),
      animate: AnimateParam(),
      distance: ReadOnlyParam(),
      radius: ReadOnlyParam(),
      sceneRadius: ReadOnlyParam(),
      showsAll: ReadOnlyParam(),
    },
    history: false,
    modules: ['canvas'],
    async run(v, ctx) {
      const canvas = ctx.plugin.canvas3d;
      if (!canvas) return ApiResult.fail('unavailable', 'There is no canvas.');
      const parsed = {
        position: parseVec(v.position),
        target: parseVec(v.target),
        up: parseVec(v.up),
      };
      for (const [name, p] of Object.entries(parsed)) {
        if (p === null) {
          return ApiResult.fail(
            'invalid',
            `'${name}' is not "x, y, z": ${JSON.stringify(v[name as keyof typeof parsed])}.`,
          );
        }
      }
      if (parsed.up && Vec3.magnitude(parsed.up) < 1e-6) return ApiResult.fail('invalid', 'The up direction is 0.');
      const current = canvas.camera.state;
      const position = parsed.position ?? current.position;
      const target = parsed.target ?? current.target;
      if (Vec3.distance(position, target) < 1e-6) {
        return ApiResult.fail(
          'invalid',
          'The camera cannot be at its target.',
          'Give a position away from the target.',
        );
      }
      const snapshot: Partial<Camera.Snapshot> = {};
      const said: string[] = [];
      if (parsed.position && !Vec3.exactEquals(parsed.position, current.position)) {
        snapshot.position = parsed.position;
        said.push(`at ${sayVec(parsed.position)}`);
      }
      if (parsed.target && !Vec3.exactEquals(parsed.target, current.target)) {
        snapshot.target = parsed.target;
        said.push(`looking at ${sayVec(parsed.target)}`);
      }
      if (parsed.up) {
        const up = Vec3.normalize(Vec3(), parsed.up);
        if (!Vec3.equals(up, current.up)) {
          snapshot.up = up;
          said.push(`up ${sayVec(up)}`);
        }
      }
      // The projection (mode, field of view) is the canvas's, set by `canvas.camera` (one operation per state).
      if (!said.length) return ApiResult.already('The camera is already so.');
      // Moving the camera: the radius it shows follows the distance to the target (as `Camera.getFocus` keeps them).
      if (snapshot.position || snapshot.target) {
        const ratio = current.radius / Math.max(Vec3.distance(current.position, current.target), 1e-6);
        snapshot.radius = Vec3.distance(position, target) * ratio;
      }
      if (Object.keys(snapshot).length) ctx.plugin.managers.camera.setSnapshot(snapshot, duration(ctx, v.animate));
      return ApiResult.ok(`Camera: ${said.join(', ')}.`, { view: true });
    },
  });

  readonly get = this.operation({
    id: 'camera.get',
    title: 'Describe the camera',
    params: {},
    read: true,
    async run(_, ctx) {
      const canvas = ctx.plugin.canvas3d;
      if (!canvas) return ApiResult.ok('No canvas (a headless viewer).', {}, { value: undefined });
      const { mode: _mode, fov: _fov, ...value } = cameraValueOf(canvas);
      const lens = cameraValueOf(canvas);
      return ApiResult.ok(
        `${lens.mode === 'perspective' ? 'Perspective' : 'Orthographic'} camera (${lens.fov}°)${value.showsAll ? ', showing the whole scene' : ', zoomed in on part of the scene'}: at ${said(value.position)}, looking at ${said(value.target)}, up ${said(value.up)}.`,
        {},
        { value },
      );
    },
  });
}

const round = (v: number) => Math.round(v * 10) / 10;
const round2 = (v: number) => Math.round(v * 100) / 100;
const vec = (v: Vec3): [number, number, number] => [round2(v[0]), round2(v[1]), round2(v[2])];
const said = (v: readonly number[]) => `(${v.join(', ')})`;
const sayVec = (v: Vec3) => said(vec(v));

/** What `camera.get` gives (and `camera.set` takes back): the view; the projection is `canvas.get`'s and `canvas.camera`'s. */
export type CameraView = Omit<CameraValue, 'mode' | 'fov'>;

/** The camera in values (the Export card's Viewport section; `camera.get` without the projection). */
export function cameraValueOf(canvas: Pick<Canvas3D, 'camera' | 'boundingSphere'>): CameraValue {
  const { mode, fov, radius, position, target, up } = canvas.camera.state;
  const sceneRadius = canvas.boundingSphere.radius;
  return {
    mode,
    fov: Math.round(radToDeg(fov)),
    position: vec(position),
    target: vec(target),
    up: vec(up),
    distance: round2(Vec3.distance(position, target)),
    radius: round(radius),
    sceneRadius: round(sceneRadius),
    showsAll: sceneRadius <= 0 || radius >= sceneRadius * 0.95,
  };
}

/** A vector `camera.set` takes: three numbers (what `camera.get` gives), or "x, y, z" text (a form's, a tool's). */
export type VecValue = string | readonly number[];

/** A vector param: a text field in a form and a string in a tool's schema; numbers from code. */
function VecParam(description: string) {
  return PD.Converted<VecValue, PD.Text>(
    (v) => (typeof v === 'string' ? v : v.join(', ')),
    (t) => t,
    PD.Text('', { description }),
  );
}

/** A field of `camera.get`'s value that `camera.set` cannot set: taken, so the value goes back as it is, and ignored. */
function ReadOnlyParam() {
  return PD.Value<unknown>(undefined, { isHidden: true, description: "camera.get's; ignored" });
}

function AnimateParam() {
  return PD.Boolean(true, { description: 'A short transition instead of a cut' });
}

/**
 * "x, y, z" (or "[x, y, z]", an array the agent gave as text), or three numbers, as a vector; undefined when empty, null
 * when not one.
 */
export function parseVec(value: unknown): Vec3 | undefined | null {
  if (Array.isArray(value)) {
    return value.length === 3 && value.every((x) => typeof x === 'number' && Number.isFinite(x))
      ? Vec3.create(value[0], value[1], value[2])
      : null;
  }
  if (value === undefined || value === null) return undefined;
  if (typeof value !== 'string') return null;
  const t = value.trim().replace(/^\[|\]$/g, '');
  if (!t) return undefined;
  const parts = t
    .split(/[\s,;]+/)
    .filter(Boolean)
    .map(Number);
  if (parts.length !== 3 || !parts.every(Number.isFinite)) return null;
  return Vec3.create(parts[0], parts[1], parts[2]);
}

/** The camera's animation, or a cut when not asked for or where the person prefers reduced motion (README §3.7). */
export function duration(ctx: { modules: { has(key: 'page'): boolean } }, animate = true): number {
  return !animate || (ctx.modules.has('page') && prefersReducedMotion()) ? 0 : ZoomDurationMs;
}

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}
