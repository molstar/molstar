/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ApiArea, ApiResult } from './area.js';
import { EveryStructureParam, StructureParam } from './params.js';
import { isFailure, labelOf, someStructures } from './scope.js';
import { LanguageParam, TargetParam } from './target-param.js';

// `camera` (operations.md §2, plan 2.3): view state, never a history step (sequence.md §7.3). `reset`, `frame` (a
// target), `zoom`, `orient` and `get`. What moves the camera needs the canvas (`modules: ['canvas']`); a headless
// viewer says so instead of running.

/** What `camera.get` gives back; undefined without a canvas. */
export interface CameraValue {
  readonly mode: 'perspective' | 'orthographic';
  readonly fov: number;
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

  readonly get = this.operation({
    id: 'camera.get',
    title: 'Describe the camera',
    params: {},
    read: true,
    async run(_, ctx) {
      const canvas = ctx.plugin.canvas3d;
      if (!canvas) return ApiResult.ok('No canvas (a headless viewer).', {}, { value: undefined });
      const { mode, fov, radius } = canvas.camera.state;
      const sceneRadius = canvas.boundingSphere.radius;
      const showsAll = sceneRadius <= 0 || radius >= sceneRadius * 0.95;
      const value: CameraValue = { mode, fov, radius: round(radius), sceneRadius: round(sceneRadius), showsAll };
      return ApiResult.ok(
        `${mode === 'perspective' ? 'Perspective' : 'Orthographic'} camera${showsAll ? ', showing the whole scene' : ', zoomed in on part of the scene'}.`,
        {},
        { value },
      );
    },
  });
}

const round = (v: number) => Math.round(v * 10) / 10;

/** Milliseconds of the transition: none when the person asks for less motion. */
/** The camera's animation, or a cut where the person prefers reduced motion (README §3.7). */
export function duration(ctx: { modules: { has(key: 'page'): boolean } }): number {
  return ctx.modules.has('page') && prefersReducedMotion() ? 0 : ZoomDurationMs;
}

function prefersReducedMotion(): boolean {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}
