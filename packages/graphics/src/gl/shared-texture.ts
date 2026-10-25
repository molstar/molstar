/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { TextureSpec } from './renderable/schema.js';
import { TextureChanges } from './renderable/texture-changes.js';
import type { TextureImage, TextureVolume } from './renderable/util.js';
import type { Texture, TextureKind } from './texture.js';

export type SharedTextureImage = TextureImage<any> | TextureVolume<any>;

/** A holder's reference to a shared texture, see `SharedTextures.acquire` */
export interface SharedTexture {
  readonly texture: Texture;
  readonly image: SharedTextureImage;
  /** Uploads the image again when its `version` changed since the last upload, always for images without one */
  sync: () => void;
  /** The holder's share of the texture's bytes (fractional): the byte counts of the items sharing it add up to it once */
  getByteCount: () => number;
  /** Destroys the texture when it was the last holder; idempotent */
  release: () => void;
}

/**
 * Textures shared by the render items of a context: one per image object and texture spec, counted by holders, for
 * specs with `shared` (`TextureSpecOptions`), see `.v6/plans/gl/shared-textures.md` §2
 */
export interface SharedTextures {
  acquire: (spec: TextureSpec<TextureKind>, image: SharedTextureImage) => SharedTexture;
  /** The number of live shared textures */
  readonly count: number;
  /**
   * Advances when an upload replaced the storage of a held texture (a grown image): draws recorded with it bind the old
   * storage, also those of the other holders, whose values did not change (`GraphicsContext.drawResourcesVersion`)
   */
  readonly version: number;
}

type Entry = {
  readonly texture: Texture;
  /** `image.version` at the last upload */
  uploaded: number | undefined;
  holders: number;
};

/** What the GPU texture depends on */
function getSpecKey(spec: TextureSpec<TextureKind>) {
  const { kind, format, dataType, filter, dimension, sampleType } = spec;
  const dim = typeof dimension === 'string' ? dimension : dimension.define;
  return `${kind}|${format}|${dataType}|${filter}|${dim}|${sampleType}`;
}

/**
 * `create` gives an unloaded texture from the backend's factory, so that resource and byte counts see it once;
 * `storage` tells a texture's storage apart (one whose bindings are recorded, e.g. WebGPU's `version`), for `version`
 */
export function createSharedTextures(
  create: (spec: TextureSpec<TextureKind>) => Texture,
  storage?: (texture: Texture) => number,
): SharedTextures {
  const entries = new WeakMap<SharedTextureImage, Map<string, Entry>>();
  let count = 0;
  let version = 0;

  /** All of `image`, or its `rows` (`Texture.loadRows`) */
  function load(entry: Entry, image: SharedTextureImage, rows?: readonly number[]) {
    const { texture } = entry;
    const before = storage?.(texture);
    if (rows) texture.loadRows!(image as TextureImage<any>, rows);
    else texture.load(image);
    entry.uploaded = (image as TextureImage<any>).version;
    // a texture without holders yet was bound by no draw
    if (entry.holders > 0 && before !== storage?.(texture)) version += 1;
  }

  /** Uploads what changed since the last upload: the rows written in place where logged (`TextureChanges.writeRows`) */
  function sync(entry: Entry, image: SharedTextureImage) {
    const current = (image as TextureImage<any>).version;
    if (current !== undefined && current === entry.uploaded) return;
    const rows =
      current !== undefined && entry.uploaded !== undefined && entry.texture.loadRows && !('depth' in image)
        ? TextureChanges.getImageRows(image, entry.uploaded, current)
        : undefined;
    load(entry, image, rows);
  }

  return {
    acquire: (spec, image) => {
      const key = getSpecKey(spec);
      let byImage = entries.get(image);
      if (!byImage) entries.set(image, (byImage = new Map()));
      let entry = byImage.get(key);
      if (entry) {
        // another holder may have uploaded an earlier state of the image
        sync(entry, image);
      } else {
        entry = { texture: create(spec), uploaded: undefined, holders: 0 };
        try {
          load(entry, image);
        } catch (e) {
          entry.texture.destroy();
          throw e;
        }
        byImage.set(key, entry);
        count += 1;
      }
      entry.holders += 1;

      const held = entry;
      let released = false;
      return {
        texture: held.texture,
        image,
        sync: () => {
          if (!released) sync(held, image);
        },
        getByteCount: () => (released ? 0 : held.texture.getByteCount() / held.holders),
        release: () => {
          if (released) return;
          released = true;
          held.holders -= 1;
          if (held.holders > 0) return;
          held.texture.destroy();
          byImage.delete(key);
          if (byImage.size === 0) entries.delete(image);
          count -= 1;
        },
      };
    },
    get count() {
      return count;
    },
    get version() {
      return version;
    },
  };
}

/**
 * The shared texture of `spec` for `image` after a version change of its cell: synced when it is the image of `held`,
 * otherwise `held` released and the texture of `image` acquired (§4)
 */
export function updateSharedTexture(
  textures: SharedTextures,
  held: SharedTexture,
  spec: TextureSpec<TextureKind>,
  image: SharedTextureImage,
): SharedTexture {
  if (image === held.image) {
    held.sync();
    return held;
  }
  const next = textures.acquire(spec, image);
  held.release();
  return next;
}
