/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's `seqv:map` renderer (seqspec:viewer#extension-points): other tracks' overviews over the whole axis as an
// image of the tracks in miniature (the host's, `DrawView.tracks`), kept while its key holds (a pan or a zoom builds
// none), the view's window over it washed and outlined in two tones. The image is this renderer's cache, so a new
// prepared state or skin makes another.
import {
  type DrawView,
  KeyedImage,
  kernelEntry,
  SkinToken,
  WINDOW_WASH,
  windowBox,
  windowOverMap,
} from "../draw/api.js";

const img = new Float64Array(8);

function draw(v: DrawView, s: { readonly tracks: readonly string[] }): void {
  const k = v.cache as KeyedImage | null;
  const host = v.tracks;
  if (!k || !host) return;
  const d = v.dpr;
  const x = Math.round(v.left * d) | 0;
  const y = Math.round((v.top + 2) * d) | 0;
  const w = (Math.round((v.left + v.navWidth) * d) | 0) - x;
  const h = (Math.round((v.top + v.height - 2) * d) | 0) - y;
  if (w < 2 || h < 2) return;
  const images = v.backend.images;
  k.next[0] = w;
  k.next[1] = h;
  k.next[2] = d;
  k.next[3] = v.navCell;
  host.stamp(k.next, 4);
  if (k.stale(images)) {
    const drawn = host.draw(images.beginLayer(w, h), s.tracks, v, w, h);
    const id = images.endLayer();
    if (!drawn) images.dispose(id);
    k.keep(images, drawn ? id : -1);
  }
  if (k.id < 0) return;
  img[2] = img[6] = w;
  img[3] = img[7] = h;
  img[4] = x;
  img[5] = y;
  v.backend.image(k.id, img, false);
  windowBox[0] = Math.round(v.navX0 * d) | 0;
  windowBox[1] = windowBox[5] = y;
  windowBox[2] = Math.round(v.navX1 * d) | 0;
  windowBox[3] = windowBox[7] = y + h;
  windowBox[4] = x;
  windowBox[6] = x + w;
  windowOverMap(v.backend, Math.max(1, Math.round(d)) | 0, SkinToken.surface, WINDOW_WASH);
}

export const MapRenderer = kernelEntry("seqv:map", 40, {
  minHeight: 4,
  miniature: false,
  draw,
  cache: { create: () => new KeyedImage(8), dispose: (c) => (c as KeyedImage).dispose() },
});
