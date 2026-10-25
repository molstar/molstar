/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { ColorNames } from '@molstar/core/util/color/names';
import { setDebugMode } from '@molstar/core/util/debug';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import { createRenderObject } from '../render-object.js';
import { createFrame } from '../renderable.js';
import { Renderer } from '../renderer.js';
import { TextSchema } from '../renderable/text.js';
import { TextureSpec } from '../renderable/schema.js';
import { createTextureImage, type TextureImage } from '../renderable/util.js';
import { Scene } from '../scene.js';
import { createSharedTextures } from '../shared-texture.js';
import { TextureChanges } from '../renderable/texture-changes.js';
import type { Texture } from '../texture.js';
import { createContext } from '../webgl/context.js';
import { createGl } from './gl.shim.js';
import { tryGetGLContext } from './gl.js';

// .v6/plans/gl/shared-textures.md

/** An atlas-like image: one object, written in place, with a `version` */
function createAtlas(n = 4096) {
  const image = createTextureImage(n, 1, Uint8Array) as TextureImage<Uint8Array>;
  image.version = 0;
  return image;
}

/** A text of one glyph quad sampling `atlas`, so that drawing it binds the texture */
function createText(atlas: TextureImage<Uint8Array>) {
  const text = Text.create(
    atlas,
    new Float32Array(12),
    new Float32Array([-1, -1, -1, 1, 1, 1, 1, -1]),
    new Float32Array(4),
    new Uint32Array([0, 1, 2, 0, 2, 3]),
    new Float32Array(4),
    new Float32Array([0, 0, 0, 1, 1, 1, 1, 0]),
    1,
  );
  const props = PD.getDefaultValues(Text.Params);
  const values = Text.Utils.createValuesSimple(text, props, ColorNames.orange, 1);
  return createRenderObject('text', values, Text.Utils.createRenderableState(props), -1);
}

/** Uploads of `image`'s size: whole (`texImage2D`, width and height at 3 and 4) or sub-image (`texSubImage2D`, at 4 and 5) */
function countUploads(gl: WebGLRenderingContext, image: TextureImage<any>) {
  let count = 0;
  const wrap = (name: 'texImage2D' | 'texSubImage2D', at: number) => {
    const f = gl[name] as (...args: any[]) => any;
    (gl as any)[name] = (...args: any[]) => {
      if (args[at] === image.width && args[at + 1] === image.height) count += 1;
      return f.apply(gl, args);
    };
  };
  wrap('texImage2D', 3);
  wrap('texSubImage2D', 4);
  return () => count;
}

describe('shared textures, the cache', () => {
  function fakeTexture() {
    return { load: vi.fn(), destroy: vi.fn(), getByteCount: () => 100 } as unknown as Texture & {
      load: ReturnType<typeof vi.fn>;
      destroy: ReturnType<typeof vi.fn>;
    };
  }
  const spec = TextSchema.tFont;

  it('one texture per image and spec, destroyed with its last holder', () => {
    const created: ReturnType<typeof fakeTexture>[] = [];
    const textures = createSharedTextures(() => {
      created.push(fakeTexture());
      return created[created.length - 1];
    });
    const atlas = createAtlas(16);
    const a = textures.acquire(spec, atlas);
    const b = textures.acquire(spec, atlas);
    expect(a.texture).toBe(b.texture);
    expect(textures.count).toBe(1);
    expect(created[0].load).toHaveBeenCalledTimes(1);
    expect(a.getByteCount() + b.getByteCount()).toBe(100);

    // another spec of the same image is another texture
    const nearest = TextureSpec('image-uint8', 'r', 'ubyte', 'nearest', undefined, { shared: true });
    const c = textures.acquire(nearest, atlas);
    expect(c.texture).not.toBe(a.texture);
    c.release();

    a.release();
    a.release(); // idempotent
    expect(created[0].destroy).not.toHaveBeenCalled();
    b.release();
    expect(created[0].destroy).toHaveBeenCalledTimes(1);
    expect(textures.count).toBe(0);
  });

  it('uploads again only when the version changed, images without one always', () => {
    const texture = fakeTexture();
    const textures = createSharedTextures(() => texture);
    const atlas = createAtlas(16);
    const a = textures.acquire(spec, atlas);
    const b = textures.acquire(spec, atlas);
    a.sync();
    b.sync();
    expect(texture.load).toHaveBeenCalledTimes(1);
    atlas.version! += 1;
    a.sync();
    b.sync();
    expect(texture.load).toHaveBeenCalledTimes(2);
    // a later holder of an image changed since the last upload
    atlas.version! += 1;
    textures.acquire(spec, atlas);
    expect(texture.load).toHaveBeenCalledTimes(3);

    const plain = createTextureImage(16, 1, Uint8Array);
    const p = textures.acquire(spec, plain);
    p.sync();
    expect(texture.load).toHaveBeenCalledTimes(5);
  });

  it('advances its version when an upload replaces the storage of a held texture', () => {
    // a texture whose storage is replaced when the image's size changes, as WebGPU's (`WebGPUTexture.version`)
    let storage = 0;
    let size = -1;
    const texture = {
      ...fakeTexture(),
      load: (image: { width: number }) => {
        if (image.width !== size) storage += 1;
        size = image.width;
      },
    } as unknown as Texture;
    const textures = createSharedTextures(
      () => texture,
      () => storage,
    );
    const atlas = createAtlas(16);
    // the first upload: nothing was bound to it
    const a = textures.acquire(spec, atlas);
    textures.acquire(spec, atlas);
    expect(textures.version).toBe(0);
    // new glyphs in place: the same storage
    atlas.version! += 1;
    a.sync();
    expect(textures.version).toBe(0);
    // the atlas grown: the other holder's recorded draws are stale too
    (atlas as { width: number }).width = 32;
    atlas.version! += 1;
    a.sync();
    expect(textures.version).toBe(1);
  });

  it('uploads the rows written since its last upload, all where they are not logged', () => {
    const texture = { ...fakeTexture(), loadRows: vi.fn() } as unknown as Texture & {
      load: ReturnType<typeof vi.fn>;
      loadRows: ReturnType<typeof vi.fn>;
    };
    const textures = createSharedTextures(() => texture);
    const atlas = createAtlas(64 * 64);
    const a = textures.acquire(spec, atlas);
    expect(texture.load).toHaveBeenCalledTimes(1);
    // two glyphs in the same rows, one below
    TextureChanges.writeRows(atlas, 0, 4);
    TextureChanges.writeRows(atlas, 0, 4);
    TextureChanges.writeRows(atlas, 4, 8);
    a.sync();
    expect(texture.loadRows).toHaveBeenCalledTimes(1);
    expect(texture.loadRows).toHaveBeenLastCalledWith(atlas, [0, 8]);
    expect(texture.load).toHaveBeenCalledTimes(1);
    // written without rows: all of it
    atlas.version! += 1;
    a.sync();
    expect(texture.load).toHaveBeenCalledTimes(2);
    expect(texture.loadRows).toHaveBeenCalledTimes(1);
  });

  it("kind 'texture' cannot be shared", () => {
    expect(() => TextureSpec('texture', 'rgba', 'ubyte', 'linear', undefined, { shared: true })).toThrow();
  });
});

describe('shared textures, WebGL render items', () => {
  function setup() {
    const gl = createGl(32, 32, { preserveDrawingBuffer: true });
    const ctx = createContext(gl);
    const scene = Scene.create(ctx);
    return { gl, ctx, scene };
  }

  it('two texts with the same atlas share one texture, uploaded once', () => {
    const { gl, ctx, scene } = setup();
    const atlas = createAtlas();
    const uploads = countUploads(gl, atlas);
    const textures = ctx.stats.resourceCounts.texture;
    const a = createText(atlas);
    const b = createText(atlas);
    scene.add(a);
    scene.add(b);
    scene.commit();
    expect(ctx.resources.sharedTextures.count).toBe(1);
    expect(uploads()).toBe(1);

    // a third text with its own atlas: its own texture
    const other = createAtlas(1024);
    const c = createText(other);
    scene.add(c);
    scene.commit();
    expect(ctx.resources.sharedTextures.count).toBe(2);

    // the scene's renderables count a shared texture once
    const atlasBytes = atlas.width * atlas.height;
    const [ra, rb] = scene.renderables.filter((r) => r.values.tFont.ref.value === atlas);
    expect(ra.getByteCount() + rb.getByteCount()).toBeGreaterThanOrEqual(atlasBytes);
    expect(ra.getByteCount() + rb.getByteCount()).toBeLessThan(2 * atlasBytes);
    // one per atlas for tFont, one per item for its other image texture (the palette), the items' data arrays in the
    // data pages (a texture per format, `data-pages.ts`)
    const pages = ctx.stats.resourceCounts.texture - textures - 2 - 3;
    expect(pages).toBeGreaterThan(0);
    expect(pages).toBeLessThanOrEqual(3);
    ctx.destroy();
  });

  it('an update uploads once for all items, and not at all without new glyphs', () => {
    const { gl, ctx, scene } = setup();
    const atlas = createAtlas();
    const a = createText(atlas);
    const b = createText(atlas);
    scene.add(a);
    scene.add(b);
    scene.commit();
    const uploads = countUploads(gl, atlas);

    // a rebuild: both cells updated with the same image, no glyph appended
    ValueCell.update(a.values.tFont, atlas);
    ValueCell.update(b.values.tFont, atlas);
    scene.update([a, b], false);
    expect(uploads()).toBe(0);

    // a glyph appended in place
    atlas.array[0] = 255;
    atlas.version! += 1;
    ValueCell.update(a.values.tFont, atlas);
    ValueCell.update(b.values.tFont, atlas);
    scene.update([a, b], false);
    expect(uploads()).toBe(1);
    expect(ctx.resources.sharedTextures.count).toBe(1);

    // a glyph written as rows (the font atlas): those rows uploaded, once for all items
    const rows: number[] = [];
    const texSubImage2D = gl.texSubImage2D.bind(gl);
    (gl as any).texSubImage2D = (...args: any[]) => {
      rows.push(args[3], args[5]);
      return (texSubImage2D as any)(...args);
    };
    atlas.array.fill(7, 2 * atlas.width, 3 * atlas.width);
    TextureChanges.writeRows(atlas, 2, 3);
    ValueCell.update(a.values.tFont, atlas);
    ValueCell.update(b.values.tFont, atlas);
    scene.update([a, b], false);
    expect(rows).toEqual([2, 1]);
    expect(uploads()).toBe(1);
    ctx.destroy();
  });

  it('a new atlas object releases the old texture and acquires its own', () => {
    const { gl, ctx, scene } = setup();
    const created = vi.spyOn(gl, 'createTexture');
    const deleted = vi.spyOn(gl, 'deleteTexture');
    const atlas = createAtlas();
    const a = createText(atlas);
    const b = createText(atlas);
    scene.add(a);
    scene.add(b);
    scene.commit();
    const old = created.mock.results.map((r) => r.value);

    const next = createAtlas(1024);
    const uploads = countUploads(gl, next);
    ValueCell.update(a.values.tFont, next);
    scene.update([a], false);
    expect(ctx.resources.sharedTextures.count).toBe(2);
    ValueCell.update(b.values.tFont, next);
    scene.update([b], false);
    expect(ctx.resources.sharedTextures.count).toBe(1);
    // one upload of the new atlas for both, and the old atlas's texture deleted (the draw after the swap is in the
    // native test below: a list still holding the deleted texture fails the bind there)
    expect(uploads()).toBe(1);
    const deletedTextures = deleted.mock.calls.map((c) => c[0]);
    expect(deletedTextures.filter((t) => old.includes(t)).length).toBe(1);
    ctx.destroy();
  });

  it('removing one keeps the texture for the other, removing both frees it', () => {
    const { gl, ctx, scene } = setup();
    const deleted = vi.spyOn(gl, 'deleteTexture');
    const atlas = createAtlas();
    const a = createText(atlas);
    const b = createText(atlas);
    scene.add(a);
    scene.add(b);
    scene.commit();
    const textures = ctx.stats.resourceCounts.texture;

    scene.remove(a);
    scene.commit();
    expect(ctx.resources.sharedTextures.count).toBe(1);
    const perItem = textures - ctx.stats.resourceCounts.texture;
    const deletedWithOne = deleted.mock.calls.length;

    // the remaining item still uploads into the shared texture
    const uploads = countUploads(gl, atlas);
    atlas.version! += 1;
    ValueCell.update(b.values.tFont, atlas);
    scene.update([b], false);
    expect(uploads()).toBe(1);

    scene.remove(b);
    scene.commit();
    expect(ctx.resources.sharedTextures.count).toBe(0);
    // the item's own textures and the atlas
    expect(deleted.mock.calls.length - deletedWithOne).toBe(perItem + 1);
    ctx.destroy();
  });
});

describe('shared textures, native WebGL', () => {
  const probe = tryGetGLContext(32, 32);
  probe?.destroy();

  (probe ? it : it.skip)(
    'a text drawn after the other holder of its atlas is gone, and after a swap',
    { tags: ['medium'] },
    () => {
      const ctx = tryGetGLContext(32, 32)!;
      const gl = ctx.gl;
      const scene = Scene.create(ctx);
      const renderer = Renderer.create(ctx);
      const camera = new Camera({ position: Vec3.create(0, 0, 50) });
      const atlas = createAtlas();
      const objects = [createText(atlas), createText(atlas)];
      for (const o of objects) scene.add(o);
      scene.commit();
      expect(ctx.resources.sharedTextures.count).toBe(1);

      scene.remove(objects[0]);
      scene.commit();
      atlas.version! += 1;
      ValueCell.update(objects[1].values.tFont, atlas);
      scene.update([objects[1]], false);

      expect(scene.renderables[0].values.drawCount.ref.value).toBe(6);
      ctx.resources.linkPrograms();
      expect(ctx.resources.finalizePrograms(undefined, true)).toBe(true);
      setDebugMode(true);
      try {
        renderer.update(camera, scene, createFrame());
        ctx.beginRenderPass('drawing-buffer');
        expect(() => renderer.renderBlended(scene, camera)).not.toThrow();
        ctx.endRenderPass();
        expect(gl.getError()).toBe(gl.NO_ERROR);

        // swapped to another atlas, the last holder of the first: its texture is deleted, and a draw that still bound it
        // would raise INVALID_OPERATION
        const next = createAtlas(1024);
        ValueCell.update(objects[1].values.tFont, next);
        scene.update([objects[1]], false);
        expect(ctx.resources.sharedTextures.count).toBe(1);
        renderer.update(camera, scene, createFrame());
        ctx.beginRenderPass('drawing-buffer');
        expect(() => renderer.renderBlended(scene, camera)).not.toThrow();
        ctx.endRenderPass();
        expect(gl.getError()).toBe(gl.NO_ERROR);
      } finally {
        setDebugMode(false);
      }

      scene.remove(objects[1]);
      scene.commit();
      expect(ctx.resources.sharedTextures.count).toBe(0);
      renderer.dispose();
      ctx.destroy();
    },
  );
});
