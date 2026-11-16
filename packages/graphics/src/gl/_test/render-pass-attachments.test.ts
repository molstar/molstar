/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { type BackendRenderPassAttachments, createRenderPassAttachments } from '../render-pass-attachments.js';
import { createWebGPUContext } from '../webgpu/context.js';
import type { WebGPUAttachmentSet } from '../webgpu/attachment-set.js';
import { createStubDevice } from './stub-gpu.js';

function createCase() {
  const ctx = createWebGPUContext(createStubDevice(), { width: 16, height: 16 });
  const color = () => {
    const t = ctx.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    t.define(16, 16);
    return t;
  };
  const depth = ctx.resources.texture('image-depth', 'depth', 'ushort', 'nearest');
  depth.define(16, 16);
  return { ctx, color, depth };
}

/** the attachment set the backend draws them with */
const getAttachmentSet = (attachments: unknown) =>
  (attachments as BackendRenderPassAttachments<WebGPUAttachmentSet>).getBackendTarget();

const refs = (fb: { getAttachments(): readonly { kind: string; id: number }[] }) =>
  fb.getAttachments().map((r) => `${r.kind}:${r.id}`);

// a pass names its targets as a whole (as WebGPU's render pass descriptor): the same again is no change, any field
// changed attaches it
describe('render pass attachments', () => {
  it('draws into what they name, a pass again changing nothing', () => {
    const { ctx, color, depth } = createCase();
    const a = color();
    const b = color();
    const attachments = createRenderPassAttachments(ctx.resources, [a, b], depth);

    ctx.beginRenderPass(attachments);
    ctx.endRenderPass();
    const set = getAttachmentSet(attachments);
    expect(refs(set)).toEqual([`texture:${a.id}`, `texture:${b.id}`, `texture:${depth.id}`]);
    const version = set.version;
    ctx.beginRenderPass(attachments);
    ctx.endRenderPass();
    expect(set.version).toBe(version);

    // a texture swapped (ping-pong), one color fewer, no depth: attached as named
    const c = color();
    attachments.color = [c];
    attachments.depth = undefined;
    ctx.beginRenderPass(attachments);
    ctx.endRenderPass();
    expect(set.version).toBeGreaterThan(version);
    expect(refs(set)).toEqual([`texture:${c.id}`]);
    ctx.destroy();
  });

  it('describes the target of its formats, and is a resource', () => {
    const { ctx, color, depth } = createCase();
    const attachments = createRenderPassAttachments(ctx.resources, [color()], depth);
    expect(ctx.describeTarget(attachments).key).toBe(ctx.describeTarget(attachments).key);
    attachments.destroy();
    expect(() => getAttachmentSet(attachments)).toThrow(/destroyed/);
    ctx.destroy();
  });
});
