/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { decodeImage, decodeImageBitmap } from '../image.js';

// The browser's image decoder: none in Node, refused there by the image's name (a reader takes a decoder of its own)

describe('decodeImage', () => {
  it('refuses outside a browser, naming the image (Node has no decoder)', async () => {
    const bytes = new Uint8Array(4);
    await expect(decodeImage(bytes, { type: 'image/webp', name: 'sh0.webp' })).rejects.toThrow(
      "decoding the image 'sh0.webp' needs a browser",
    );
    expect(() => decodeImageBitmap(bytes, 'image/png', 'a.png')).toThrow("decoding the image 'a.png' needs a browser");
  });
});
