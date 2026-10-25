/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type {
  GraphicsContext,
  PipelineRequest,
  RenderPassOptions,
  RenderPassTarget,
  TargetDescriptor,
} from '@molstar/graphics/gl/context';
import type { ComputeRenderable } from '@molstar/graphics/gl/renderable';
import type { RenderState } from '@molstar/graphics/gl/render-state';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Scene } from '@molstar/graphics/gl/scene';

/**
 * What the passes draw, written once (their `plan` methods) and interpreted twice: rendering draws it, preparing
 * collects the pipelines it draws with (`resources.preparePipelines`) without drawing anything.
 *
 * A plan runs its control flow under both. Code that only matters for drawing (per-frame uniforms, viewports, clears,
 * timers, readbacks) goes behind `drawing`; code that selects programs or targets (props applied to renderables, e.g.
 * defines, and attachment changes, which change the formats of targets) runs in both. Every target a frame draws into
 * is begun in both, also for a clear only: plans name the targets that keep their storage (`PassTargets`).
 */
export interface PassBuilder {
  /** true when rendering, false when preparing */
  readonly drawing: boolean;
  /** starts a render pass on `target` (preparing: the target the following draws are collected for) */
  begin(target: RenderPassTarget, options?: RenderPassOptions): void;
  /**
   * a renderer pass over `group`, `draw` draws it (`collect` describes it); `all`: prepare also the objects the pass
   * leaves out in their current state (e.g. unmarked objects for the marking mask)
   */
  renderer(pass: Renderer.Pass, group: Scene.Group, draw: () => void, all?: boolean): void;
  /** a full-screen renderable drawn with `state` */
  quad(renderable: ComputeRenderable<any>, state: RenderState): void;
}

/** Draws the plan */
export function createRenderBuilder(graphics: GraphicsContext): PassBuilder {
  return {
    drawing: true,
    begin: (target, options) => graphics.beginRenderPass(target, options),
    renderer: (_pass, _group, draw) => draw(),
    quad: (renderable, state) => {
      graphics.setRenderState(state);
      renderable.render();
    },
  };
}

/**
 * Requests the pipelines of the plan as it runs: a render item's program depends on its current defines and textures,
 * which a plan may change between draws (e.g. a blur kernel size per mip level), so each request is made right away
 */
export interface PrepareBuilder extends PassBuilder {
  /** the requests made so far */
  readonly requests: PipelineRequest[];
  /** resolves when the requested pipelines (and other pending work) are ready */
  whenReady(): Promise<void>;
}

/**
 * The requests of a renderer pass over `group`, `undefined` for the renderer's own (`renderer.collect`); e.g. for
 * objects described by signatures instead of the group's renderables (`Canvas3D.precompile`)
 */
export type PrepareCollect = (
  pass: Renderer.Pass,
  group: Scene.Group,
  target: TargetDescriptor,
  all: boolean,
) => PipelineRequest[] | undefined;

export function createPrepareBuilder(
  graphics: GraphicsContext,
  renderer: Renderer,
  collect?: PrepareCollect,
): PrepareBuilder {
  let target: TargetDescriptor | undefined;
  const getTarget = () => {
    if (!target) throw new Error('draw outside of a render pass');
    return target;
  };
  const requests: PipelineRequest[] = [];
  const request = (r: PipelineRequest[]) => {
    for (const x of r) requests.push(x);
    graphics.resources.preparePipelines(r);
  };
  return {
    drawing: false,
    requests,
    whenReady: () => graphics.resources.pending.whenReady(),
    begin: (t, options) => {
      target = graphics.describeTarget(t, options);
    },
    renderer: (pass, group, _draw, all) => {
      const target = getTarget();
      request(collect?.(pass, group, target, !!all) ?? renderer.collect(pass, group, target, { all }));
    },
    quad: (renderable, state) =>
      request([{ item: renderable.renderItem, variant: 'compute', state, target: getTarget() }]),
  };
}
