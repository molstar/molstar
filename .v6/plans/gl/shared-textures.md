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

- Uploads are of the whole atlas: a new glyph uploads 43 MB once. A sub-image upload of the rows written since the last
  upload would cut that; the atlas would record a dirty range next to `version`.
- WebGPU render items are not covered by tests (no WebGPU context in the test lanes): acquire, the swap that resolves
  bind groups again, release and the byte share were checked by hand in Chrome only.
- A render item whose constructor throws after `createTextures` keeps its holders: the texture is never freed before
  its context (per-item textures leaked the same way before; a shared one pins the atlas).
- The native test checks a draw raises no GL error, not the pixels: a skipped upload after a version bump would pass.
