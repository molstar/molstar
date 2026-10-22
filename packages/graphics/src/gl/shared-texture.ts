/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TextureSpec } from './renderable/schema.js';
import type { TextureImage, TextureVolume } from './renderable/util.js';
import type { Texture, TextureKind } from './texture.js';

export type SharedTextureImage = TextureImage<any> | TextureVolume<any>;

/** A holder's reference to a shared texture, see `SharedTextures.acquire` */
export interface SharedTexture {
  readonly texture: Texture;
  readonly image: SharedTextureImage;
  /** Uploads the image again when its `version` changed since the last upload, always for images without one */
  sync: () => void;
  /** The holder's share of the texture's bytes: the byte counts of the items sharing it add up to it once */
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
}

type Entry = {
  readonly texture: Texture;
  /** `image.version` at the last upload */
  uploaded: number | undefined;
  holders: number;
};

function getSpecKey(spec: TextureSpec<TextureKind>) {
  return `${spec.kind}|${spec.format}|${spec.dataType}|${spec.filter}`;
}

/** `create` gives an unloaded texture from the backend's factory, so that resource and byte counts see it once */
export function createSharedTextures(create: (spec: TextureSpec<TextureKind>) => Texture): SharedTextures {
  const entries = new WeakMap<SharedTextureImage, Map<string, Entry>>();
  let count = 0;

  function load(entry: Entry, image: SharedTextureImage) {
    entry.texture.load(image);
    entry.uploaded = (image as TextureImage<any>).version;
  }

  function sync(entry: Entry, image: SharedTextureImage) {
    const version = (image as TextureImage<any>).version;
    if (version === undefined || version !== entry.uploaded) load(entry, image);
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
        load(entry, image);
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
