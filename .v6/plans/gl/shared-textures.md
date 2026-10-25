# Plan: Shared image textures (one GPU copy of the font atlas)

## Status

- **Current**: done for `tFont` on WebGL and WebGPU. Other image textures keep a texture per render item.
- Why: each render item created and uploaded a texture of its own per image texture (`createTextures` of both
  backends). Every text render object (a measurement, an MVS labelled structure or primitives node, a label
  representation, the camera helper's axes) held its own copy of the font atlas, 43 MB at `fontQuality` 3, and uploaded
  all of it again on every rebuild (`Text.create` bumps `tFont`'s cell). The atlas is one per font
  ([labels.md](../../designs/labels.md) §1), so the copies are identical.

## Design

1. **Opt-in per texture spec**: `TextureSpecOptions.shared` (`renderable/schema.ts`), set for `tFont`
   (`renderable/text.ts`). Kinds other than `'texture'` only (those are owned outside the render item). Why opt-in:
   most image textures (`tColor`, `tMarker`, `tSize`, …) are per object and written per item, so sharing buys nothing
   and a lookup costs; specs whose images are shared objects opt in.
2. **Cache**: `createSharedTextures` (`gl/shared-texture.ts`, backend-neutral), one per context: `resources.sharedTextures`
   on WebGL, `WebGPUDrawContext.sharedTextures` on WebGPU. Keyed by the image object (`TextureImage`,
   `TextureVolume`) in a `WeakMap`, then by the spec's kind, format, data type and filter (the GPU texture depends on
   them). `acquire` creates and loads on first use and counts holders; `release` decrements and destroys the texture at
   zero. The texture comes from the backend's usual factory (`resources.texture`, `createWebGPUTexture`), so resource
   counts and `getByteCounts` see it once.
3. **Content version**: `TextureImage.version`, incremented by whoever writes `array` in place (`FontAtlas.get`
   appending a glyph). `sync` uploads only when it differs from the version of the last upload, so N items whose cells
   changed upload once, and a rebuild with no new glyph uploads nothing, when it updates the render objects in place (a
   visual that replaces its only text object releases the texture to zero holders first, and the new one uploads it
   again). Images without a version upload on every
   `acquire` and `sync`: correct, shared memory, no saved uploads.
   - **Rows**: a writer that knows the rows it wrote logs them with the version (`TextureChanges.writeRows`; the atlas
     its glyph's strip), and `sync` uploads the rows written since its last upload (`TextureChanges.getImageRows`,
     merged; `Texture.loadRows`, a sub-image write per range on both backends) instead of the image. A version without
     rows, a holder further behind than the log, or rows covering more than half the image upload all of it. A glyph
     uploads its strip (a line height of rows), not the 43 MB atlas at `fontQuality` 3.
   - **Fixed size**: the atlas is allocated once at its size and gives the placeholder glyph when full, so uploads never
     replace its storage (§7).
4. **Render items** (`webgl/render-item.ts`, `webgpu/render-item.ts`): `createTextures` acquires the shared ones. On a
   version change of the cell: the same image object → `sync`; another object → release the old, acquire the new and
   replace it in the item's texture list (WebGPU: as a replaced texture, the bind groups and linear-read programs are
   resolved again). `destroy` releases instead of destroying.
5. **Byte counts**: `RenderItem.getByteCount` counts a shared texture's bytes divided by its holders, so the sum over
   all its holders counts it once. Holders are per context, so a sum over one scene (Canvas3D's main scene, without the
   helper scene's camera axes) can come out a share short; `resources.getByteCounts` sees one texture.
6. **Context loss**: restore is a rebuild on a new context (`GraphicsContext.recreate`,
   [backend-switch-followups.md](backend-switch-followups.md)): the new context has a new cache, the render items of the old one are destroyed and
   release into the old cache, whose textures are gone with their context (the WebGL wrapper skips deleting them on a
   lost context).
7. **Draw lists** (WebGPU render bundles): a recording binds the texture's storage. An upload that replaces it (an image
   of another size) leaves the other holders' recordings bound to the old one, while their values did not change. So
   the cache counts storage replacements of held textures (`SharedTextures.version`, from the backend's storage
   identity: `WebGPUTexture.version`), and `GraphicsContext.drawResourcesVersion` adds it: recordings made before do
   not replay, also those of chunks whose rows did not change (the chunk check compares it before its view fallback,
   `isChunkInputsCurrent`). The font atlas keeps its size (a full atlas gives the placeholder glyph), so labels never replace it.

## Measurements

Prism viewer (`prism-viewer`), 1280 × 1536 drawing buffer, `resources.getByteCounts().texture` after the scene
settled, before (`shared` off on `tFont`) → after; the same on WebGL and WebGPU:

| Scene                                                     | Text objects | Before   | After    |
| --------------------------------------------------------- | ------------ | -------- | -------- |
| 1cbs, 10 distance measurements                            | 10           | 983.7 MiB | 613.8 MiB |
| MVS story, 6 structures with 3 labels each (1st snapshot) | 6            | 819.3 MiB | 613.8 MiB |

The saving is (text objects − 1) × 41.1 MiB, one atlas at `fontQuality` 3. A rebuild of the texts in place without new
glyphs uploads nothing (it uploaded the atlas once per text object).

## Open

- WebGPU render items are covered by native tests (`plugin/headless/src/_test/webgpu-shared-textures.test.ts`, the
  Dawn lane): two label texts sharing the atlas, new glyphs uploaded as rows, the GPU atlas read back and compared,
  storage and `drawResourcesVersion` kept; 65 texts (two chunks of a draw list) with the atlas grown under one of them,
  recorded again with no WebGPU error; a text given another atlas (its texture acquired, the old one freed with its
  last holder) and the byte shares. `Texture.loadRows` is checked on every profile (`gpu-paths/load-rows`, visual
  suite).
- A render item whose constructor throws after `createTextures` keeps its holders: the texture is never freed before
  its context (per-item textures leaked the same way before; a shared one pins the atlas).
- The native test checks a draw raises no GL error, not the pixels: a skipped upload after a version bump would pass.
