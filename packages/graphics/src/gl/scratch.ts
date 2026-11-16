/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from './context.js';

/**
 * Scratch resources of one module per context (textures, framebuffers, compute renderables a computation fills and
 * reads, reused by the next): by name within the module, so that modules can't collide; released with the context's
 * resources when it is destroyed (`GraphicsResources.destroy`), not before
 */
export interface Scratch<T> {
  /** what `name` holds for `ctx`, undefined before it was set */
  get(ctx: GraphicsContext, name?: string): T | undefined;
  /** `value` held for `ctx` under `name`, returned */
  set(ctx: GraphicsContext, value: T, name?: string): T;
  /** what `name` holds for `ctx`, made by `create` the first time */
  getOrCreate(ctx: GraphicsContext, create: () => T, name?: string): T;
}

/** The scratch resources of a module; `label` names it for debugging */
export function createScratch<T>(label: string): Scratch<T> {
  const byContext = new WeakMap<GraphicsContext, Map<string, T>>();
  const getMap = (ctx: GraphicsContext) => {
    let map = byContext.get(ctx);
    if (!map) byContext.set(ctx, (map = new Map()));
    return map;
  };
  return {
    get: (ctx, name = label) => byContext.get(ctx)?.get(name),
    set: (ctx, value, name = label) => {
      getMap(ctx).set(name, value);
      return value;
    },
    getOrCreate: (ctx, create, name = label) => {
      const map = getMap(ctx);
      let value = map.get(name);
      if (value === undefined) map.set(name, (value = create()));
      return value;
    },
  };
}
