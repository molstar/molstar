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
import { EveryStructureParam, StructureParam } from './params.js';
import { isFailure, labelOf, someStructures } from './scope.js';
import { LanguageParam, TargetParam } from './target-param.js';

// `camera` (operations.md §2, plan 2.3): view state, never a history step (sequence.md §7.3). `reset`, `frame` (a
// target), `zoom`, `orient`, `set` (position, target, up, fov, mode: v5's camera info, parity.md) and `get`. What moves the camera needs the canvas (`modules: ['canvas']`); a headless
// viewer says so instead of running.

/** What `camera.get` gives back; undefined without a canvas. */
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
    params: {
      animate: PD.Boolean(true, { description: 'A short transition instead of a cut' }),
    },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ animate }, ctx) {
      const reduced = ctx.modules.has('page') && prefersReducedMotion();
      ctx.plugin.managers.camera.reset(undefined, animate && !reduced ? 250 : 0);
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
    },
    history: false,
    modules: ['canvas'],
    async run({ target, structure, language }, ctx) {
      const resolved = await ctx.targets.resolve(target, {
        ...(structure ? { structure } : {}),
        ...(language ? { language } : {}),
      });
      if (isFailure(resolved)) return resolved;
      if (resolved.empty) return ApiResult.fail('not-found', resolved.summary, 'Nothing was changed.');
      const loci = resolved.parts.filter((p) => p.selected.atoms > 0).map((p) => p.loci);
      ctx.plugin.managers.camera.focusLoci(loci.length === 1 ? loci[0] : loci, { durationMs: duration(ctx) });
      return ApiResult.ok(`Framed ${resolved.summary}.`, { view: true });
    },
  });

  readonly zoom = this.operation({
    id: 'camera.zoom',
    title: 'Zoom',
    params: {
      factor: PD.Numeric(1.5, { min: 0.1, max: 10, step: 0.1 }, { description: 'Above 1 zooms in, below 1 out' }),
    },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ factor }, ctx) {
      const canvas = ctx.plugin.canvas3d;
      if (!canvas) return ApiResult.fail('unavailable', 'There is no canvas.');
      if (factor === 1) return ApiResult.already('The zoom is unchanged.');
      // The same view direction and target, the radius it shows scaled.
      const { target, radius } = canvas.camera.state;
      const snapshot = canvas.camera.getFocus(target, Math.max(radius, 0.01) / factor);
      ctx.plugin.managers.camera.setSnapshot(snapshot, duration(ctx));
      return ApiResult.ok(factor > 1 ? `Zoomed in ${factor}×.` : `Zoomed out ${round(1 / factor)}×.`, { view: true });
    },
  });

  readonly orient = this.operation({
    id: 'camera.orient',
    title: 'Orient along the principal axes',
    params: { structure: EveryStructureParam() },
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run({ structure }, ctx) {
      const named = structure ? someStructures(ctx, structure) : undefined;
      if (named && isFailure(named)) return named;
      const data = named?.flatMap((s) => (s.cell.obj ? [s.cell.obj.data] : []));
      ctx.plugin.managers.camera.orientAxes(data, duration(ctx));
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
    params: {},
    history: false,
    modules: ['canvas'],
    when: (ctx) => (ctx.sources.empty.value ? 'Nothing is loaded.' : true),
    async run(_, ctx) {
      ctx.plugin.managers.camera.resetAxes(duration(ctx));
      return ApiResult.ok('Aligned the view with the x, y and z axes.', { view: true });
    },
  });

  // v5's camera info, settable: what `camera.get` gives, back (MolViewSpec's camera is position, target and up).
  readonly set = this.operation({
    id: 'camera.set',
    title: 'Set the camera',
    params: {
      position: PD.Text('', { description: 'Where the camera is: "x, y, z" in Å (empty: keep)' }),
      target: PD.Text('', { description: 'What it looks at: "x, y, z" in Å (empty: keep)' }),
      up: PD.Text('', { description: 'Its up direction: "x, y, z" (empty: keep)' }),
      fov: PD.Numeric(0, { min: 0, max: 130, step: 1 }, { description: 'Field of view in degrees (0: keep)' }),
      mode: PD.Select<'keep' | 'perspective' | 'orthographic'>(
        'keep',
        PD.arrayToOptions(['keep', 'perspective', 'orthographic'] as const),
      ),
      animate: PD.Boolean(true, { description: 'A short transition instead of a cut' }),
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
        if (p === null) return ApiResult.fail('invalid', `'${v[name as keyof typeof parsed]}' is not "x, y, z".`);
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
      // Mode and field of view are the canvas's props (the Canvas card follows them), as `canvas.camera` sets them.
      const lens: { mode?: Camera.Mode; fov?: number } = {};
      if (v.mode !== 'keep' && v.mode !== current.mode) {
        lens.mode = v.mode;
        said.push(v.mode);
      }
      if (v.fov > 0) {
        const fov = Math.max(10, v.fov);
        if (fov !== Math.round(radToDeg(current.fov))) {
          lens.fov = fov;
          said.push(`${fov}°`);
        }
      }
      if (!said.length) return ApiResult.already('The camera is already so.');
      // Moving the camera: the radius it shows follows the distance to the target (as `Camera.getFocus` keeps them).
      if (snapshot.position || snapshot.target) {
        const ratio = current.radius / Math.max(Vec3.distance(current.position, current.target), 1e-6);
        snapshot.radius = Vec3.distance(position, target) * ratio;
      }
      if (lens.mode || lens.fov) ctx.plugin.setCanvas3DProps({ camera: lens });
      if (Object.keys(snapshot).length) ctx.plugin.managers.camera.setSnapshot(snapshot, v.animate ? duration(ctx) : 0);
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
      const value = cameraValueOf(canvas);
      return ApiResult.ok(
        `${value.mode === 'perspective' ? 'Perspective' : 'Orthographic'} camera (${value.fov}°)${value.showsAll ? ', showing the whole scene' : ', zoomed in on part of the scene'}: at ${said(value.position)}, looking at ${said(value.target)}, up ${said(value.up)}.`,
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

/** The camera in values (`camera.get`, the Export card's Viewport section). */
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

/** "x, y, z" (or "[x, y, z]", an array the agent gave as text) as a vector; undefined when empty, null when not one. */
export function parseVec(text: string): Vec3 | undefined | null {
  const t = text.trim().replace(/^\[|\]$/g, '');
  if (!t) return undefined;
  const parts = t
    .split(/[\s,;]+/)
    .filter(Boolean)
    .map(Number);
  if (parts.length !== 3 || !parts.every(Number.isFinite)) return null;
  return Vec3.create(parts[0], parts[1], parts[2]);
}

/** Milliseconds of the transition: none when the person asks for less motion. */
/** The camera's animation, or a cut where the person prefers reduced motion (README §3.7). */
export function duration(ctx: { modules: { has(key: 'page'): boolean } }): number {
  return ctx.modules.has('page') && prefersReducedMotion() ? 0 : ZoomDurationMs;
}

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}
