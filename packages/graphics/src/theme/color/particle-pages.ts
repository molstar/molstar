/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import { Particle, type ParticleList, ParticleValues } from '@molstar/model/model/particles/particle-list';
import { ParticleSlots } from '@molstar/model/model/particles/particle-slots';
import type { ColorTheme, LocationColor } from '../color.js';

// The particle themes' colors of a slot list (`ParticleSlots`): a location's slot found in its page (the page found
// last tried first: slots in order find it without a search), its color by a painter made once a page list (each
// theme's own code over a flat list), what is shared across pages (palettes, maps) the theme's. A page's colors all
// at once (`ColorTheme.fillParticleColors`: a tile's block), a run at a time where the theme's index is runs

/** A flat list's colors: a particle's, and all of them as RGBA from `offset` */
export interface ParticlePainter {
  color(index: number): Color;
  fill(out: Uint8Array, offset: number): void;
}

/** The colors of a slot list's locations and of its pages: `painterOf` a page's list, made once a list */
export function getPagedPainting(painterOf: (page: ParticleList) => ParticlePainter, fallback: Color) {
  const painters = new Map<ParticleList, ParticlePainter>();
  const painter = (list: ParticleList) => {
    let p = painters.get(list);
    if (!p) painters.set(list, (p = painterOf(list)));
    return p;
  };
  const color: LocationColor = (location: Location): Color => {
    if (!Particle.isLocation(location)) return fallback;
    const { particles, index } = location;
    if (!particles.pages) return painter(particles).color(index);
    const p = ParticleSlots.findPage(particles, index);
    if (p < 0) return fallback;
    const page = particles.pages.pages[p];
    return painter(page.list).color(index - page.start);
  };
  const fill = (list: ParticleList, out: Uint8Array, offset: number) => painter(list).fill(out, offset);
  return { color, fill };
}

/** `f` of each page's list of a slot list, or of a flat list itself */
export function forEachPageList(particles: ParticleList, f: (list: ParticleList) => void) {
  ParticleSlots.forEachPage(particles, (page) => f(page.list));
}

/** Writes `color` as opaque RGBA at `offset` */
function setRgba(out: Uint8Array, offset: number, color: Color) {
  out[offset] = (color >> 16) & 255;
  out[offset + 1] = (color >> 8) & 255;
  out[offset + 2] = color & 255;
  out[offset + 3] = 255;
}

/**
 * The colors of `count` particles by `colorIndex` (a color of `colors` a particle, -1 `fallback`) as RGBA from
 * `offset`: runs a run at a time (a value once, its bytes copied), an array a particle at a time
 */
export function fillByIndex(
  colorIndex: ParticleValues<Int32Array> | undefined,
  colors: readonly Color[],
  count: number,
  fallback: Color,
  out: Uint8Array,
  offset: number,
) {
  if (!colorIndex) {
    for (let i = 0; i < count; ++i) setRgba(out, offset + i * 4, fallback);
    return;
  }
  if (ParticleValues.isRuns(colorIndex)) {
    const { offsets, values } = colorIndex;
    for (let r = 0, rl = values.length; r < rl; ++r) {
      const start = Math.min(offsets[r], count);
      const end = Math.min(offsets[r + 1], count);
      if (end <= start) continue;
      const ci = values[r];
      const o = offset + start * 4;
      setRgba(out, o, ci >= 0 ? colors[ci] : fallback);
      // the run's first texel repeated (doubling: few copies a run)
      let filled = 1;
      const n = end - start;
      while (filled < n) {
        const k = Math.min(filled, n - filled);
        out.copyWithin(o + filled * 4, o, o + k * 4);
        filled += k;
      }
    }
    return;
  }
  // a color's bytes once a color index: a lookup a particle
  const rgba = new Uint32Array(colors.length + 1);
  const bytes = new Uint8Array(rgba.buffer);
  setRgba(bytes, 0, fallback);
  for (let k = 0; k < colors.length; ++k) setRgba(bytes, (k + 1) * 4, colors[k]);
  if (offset % 4 === 0 && out.byteOffset % 4 === 0) {
    const words = new Uint32Array(out.buffer, out.byteOffset + offset, count);
    for (let i = 0; i < count; ++i) words[i] = rgba[colorIndex[i] + 1];
  } else {
    for (let i = 0; i < count; ++i)
      out.set(bytes.subarray((colorIndex[i] + 1) * 4, (colorIndex[i] + 2) * 4), offset + i * 4);
  }
}

/**
 * `theme`'s colors of the flat `list`'s particles as RGBA from `offset`: its own fill where it has one, else a color a
 * particle by location (a theme without it)
 */
export function fillParticleColors(theme: ColorTheme<any>, list: ParticleList, out: Uint8Array, offset: number) {
  if (theme.fillParticleColors) theme.fillParticleColors(list, out, offset);
  else fillByLocation(theme.color, list, out, offset);
}

/** `color`'s colors of the flat `list`'s particles, a location each, as RGBA from `offset` */
export function fillByLocation(color: LocationColor, list: ParticleList, out: Uint8Array, offset: number) {
  const location = Particle.Location(list, 0);
  for (let i = 0; i < list.count; ++i) {
    location.index = i;
    setRgba(out, offset + i * 4, color(location, false));
  }
}

export { setRgba as setParticleRgba };
export const DefaultParticleColor = Color(0xcccccc);
