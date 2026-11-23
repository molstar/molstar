/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Camera } from '../../camera.js';
import type { Helper } from '../../helper/helper.js';
import { OverlayHelper } from '../../helper/overlay-helper.js';
import { createFrame } from '@molstar/graphics/gl/renderable';
import { OverlayPass } from '../overlay.js';
import type { PassBuilder } from '../pass-builder.js';
import type { PassTargets } from '../pass-targets.js';
import { isRenderPassAttachments, type RenderPassAttachments } from '@molstar/graphics/gl/render-pass-attachments';

/** What the pass does with the renderer, in order */
function setup(count = 1) {
  const calls: string[] = [];
  const occlusionTest = (_: Sphere3D) => true;
  let current: typeof occlusionTest | null = occlusionTest;
  const occluders: (Renderer.Occluder | null)[] = [];
  const renderer = {
    clearDepth: () => calls.push('clearDepth'),
    update: (camera: { near: number }) => calls.push(`update near ${camera.near.toFixed(1)}`),
    prepareCullView: () => calls.push('prepareCullView'),
    renderOverlay: (_: Scene, __: Camera, occluder: Renderer.Occluder | null) => {
      calls.push(`renderOverlay occlusion ${current === null ? 'off' : 'on'}`);
      occluders.push(occluder);
    },
    setOcclusionTest: (f: typeof occlusionTest | null) => {
      current = f;
    },
    get occlusionTest() {
      return current;
    },
  } as unknown as Renderer;

  const scene = {
    count,
    renderables: [
      { state: { visible: true }, values: { boundingSphere: { ref: { value: Sphere3D.create(Vec3(), 10) } } } },
    ],
  } as unknown as Scene;
  const overlay = Object.create(OverlayHelper.prototype, {
    scene: { value: scene },
    cameras: { value: new WeakMap() },
    graphics: { value: { pixelRatio: 1 } },
    occluderDepth: { value: undefined, writable: true },
  }) as OverlayHelper;

  const depthRenderbuffer = { setSize: vi.fn(), getByteCount: () => 0 };
  const targets = {
    graphics: {
      resources: {
        renderPassAttachments: (color: unknown, depth: unknown) => ({
          kind: 'render-pass-attachments',
          color,
          depth,
          destroy: vi.fn(),
        }),
      },
    },
    renderbuffer: () => depthRenderbuffer,
  } as unknown as PassTargets;
  const texture = { id: 3 };
  const target = { texture } as unknown as RenderTarget;
  const sceneDepth = { id: 2 } as unknown as Texture;

  const begun: unknown[] = [];
  const b = {
    drawing: true,
    begin: (t: unknown) => begun.push(t),
    renderer: vi.fn(),
    quad: vi.fn(),
  } as PassBuilder;

  const camera = new Camera(
    { position: Vec3.create(0, 0, 100), radius: 5, radiusMax: 5, clipFar: true, minNear: 1 },
    { x: 0, y: 0, width: 100, height: 100 },
  );
  camera.update();
  const ctx = {
    renderer,
    camera,
    scene,
    helper: { overlay, screen: { isEnabled: false } } as unknown as Helper,
    frame: createFrame(),
  };
  return {
    pass: new OverlayPass(targets, 100, 100),
    b,
    ctx,
    calls,
    occluders,
    begun,
    texture,
    target,
    depthRenderbuffer,
    sceneDepth,
    overlay,
  };
}

describe('overlay pass', () => {
  it("draws onto the target with a depth of its own, the overlay's planes, and no occlusion test", () => {
    const { pass, b, ctx, calls, begun, texture, target, depthRenderbuffer, sceneDepth } = setup();
    pass.render(b, ctx, target, sceneDepth);
    // the target's color with the pass's own depth, named as a whole (the backends attach them when the pass begins)
    expect(begun).toHaveLength(1);
    const attachments = begun[0] as RenderPassAttachments;
    expect(isRenderPassAttachments(attachments)).toBe(true);
    expect(attachments.color).toEqual([texture]);
    expect(attachments.depth).toBe(depthRenderbuffer);
    // its objects culled on the GPU for the overlay's camera before the pass, then drawn: the overlay's bounds, 10
    // around the origin, 100 in front of the camera
    expect(calls).toEqual([
      'update near 89.9',
      'prepareCullView',
      'clearDepth',
      'update near 89.9',
      'renderOverlay occlusion off',
    ]);
    // the scene's occlusion test is back
    expect(ctx.renderer.occlusionTest).not.toBeNull();

    // the same target: the same attachments, their color not made again
    const color = attachments.color;
    pass.render(b, ctx, target, sceneDepth);
    expect(begun[1]).toBe(attachments);
    expect(attachments.color).toBe(color);
  });

  it("tests against the scene's depth with the scene's planes, and keeps it for picking", () => {
    const { pass, b, ctx, occluders, target, sceneDepth, overlay } = setup();
    pass.render(b, ctx, target, sceneDepth);
    expect(occluders).toEqual([{ depth: sceneDepth, near: ctx.camera.near, far: ctx.camera.far }]);
    expect(overlay.occluderDepth).toBe(sceneDepth);
  });

  it('does nothing without overlay objects', () => {
    const { pass, b, ctx, calls, begun, target, sceneDepth } = setup(0);
    pass.render(b, ctx, target, sceneDepth);
    expect(begun).toEqual([]);
    expect(calls).toEqual([]);
  });

  it('plans the opaque and the transparent pass when preparing (occluded parts of opaque objects)', () => {
    const { pass, b, ctx, calls, begun, target, sceneDepth } = setup();
    pass.render({ ...b, drawing: false }, ctx, target, sceneDepth);
    expect(begun).toHaveLength(1);
    expect(isRenderPassAttachments(begun[0])).toBe(true);
    expect(calls).toEqual([]);
    expect(vi.mocked(b.renderer).mock.calls.map((c) => c[0])).toEqual([
      { name: 'opaque' },
      { name: 'blendedTransparent' },
    ]);
  });
});
