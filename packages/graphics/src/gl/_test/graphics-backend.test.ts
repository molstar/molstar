/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, beforeEach, describe, expect, it, vi, type MockInstance } from 'vitest';
import { createGraphicsContext, type GraphicsBackendDescriptor } from '../backend.js';
import type { GraphicsContext } from '../context.js';

const canvas = {} as HTMLCanvasElement;
const context = { capabilities: { backend: 'webgl2' } } as unknown as GraphicsContext;

function backend(id: string, create: () => Promise<GraphicsContext | undefined>): GraphicsBackendDescriptor {
  return { id, isAvailable: () => true, create } as unknown as GraphicsBackendDescriptor;
}

describe('createGraphicsContext', () => {
  let warn: MockInstance;
  beforeEach(() => {
    warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
  });
  afterEach(() => warn.mockRestore());

  it('tries the next backend when one throws (e.g. a failed device request)', async () => {
    const failing = backend('webgpu', async () => {
      throw new Error('device request failed');
    });
    const working = backend('webgl2', async () => context);
    const created = await createGraphicsContext(canvas, [failing, working], {});
    expect(created?.backend.id).toBe('webgl2');
    expect(created?.context).toBe(context);
    expect(warn).toHaveBeenCalledTimes(1);
  });

  it('tries the next backend when one creates no context', async () => {
    const unavailable = backend('webgpu', async () => undefined);
    const working = backend('webgl2', async () => context);
    expect((await createGraphicsContext(canvas, [unavailable, working], {}))?.backend.id).toBe('webgl2');
  });

  it('gives undefined when no backend creates a context', async () => {
    const failing = backend('webgpu', async () => {
      throw new Error('device request failed');
    });
    expect(await createGraphicsContext(canvas, [failing], {})).toBeUndefined();
  });
});
