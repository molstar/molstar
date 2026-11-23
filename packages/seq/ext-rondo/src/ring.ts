/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type Address,
  type RefValue,
  type SectionBuilder,
  theme,
  type ViewBuilder,
} from '@molstar/seqx-core/extension';

// `rondoRing` (_spec/README.md, The builder helper): the genome ring spelled, never inferred: the ring section of an
// assembly (`center: none`: the ribbons own the hole), its cytobands and the hue ring (two `bands` layers over the same
// data), the `rondo:ribbons` layer with `fallback: 'pairs'` and its colour written (so a reader without the extension
// draws the same bins in the same colours), and the square navigator over the ring (`of: [ring, ring]`,
// `switch: include`). The spec's `extensions` is the writer's: `createSeqSpec({ extensions: [rondo] })`.

export interface RondoRingOptions {
  /** The section's id (`ring`). */
  readonly id?: string;
  readonly assembly: RefValue<'assembly'>;
  readonly links: RefValue<'pairs'>;
  /** The cytobands (a BED over the assembly with a `stain` field): the bands and the hue ring. */
  readonly bands?: RefValue<'ranges'>;
  /** A second ring under the first, its `include` the one a click opens (`below`); one ring (`inplace`, the default). */
  readonly detail?: 'inplace' | 'below';
  /** The square navigator over the ring (true by default). */
  readonly square?: boolean;
  /** The 24 hues (`rondo:hues`, not colour-blind safe) instead of the safe `categories` (false by default). */
  readonly hues?: boolean;
  /** A double-click on a ribbon, stem or stub opens it in the detail ring (or in place): `opens` written. */
  readonly opens?: boolean;
  /** The ring's diameter, px; the width's when absent. */
  readonly height?: number;
}

/** What `rondoRing` wrote: the ring section, the detail ring (`below`), the ribbons' track. */
export interface RondoRing {
  readonly ring: SectionBuilder;
  readonly detail?: SectionBuilder;
}

/** The ribbons' colour: a hue per chromosome, its `a` end's system, in the assembly's order (taken before `include`). */
export function ribbonsColor(hues: boolean) {
  return hues
    ? theme.scheme('rondo:hues', { field: 'system' })
    : theme.field('system', { palette: theme.categorical({ colors: 'categories', repeat: true }) });
}

function ringTracks(s: SectionBuilder, o: RondoRingOptions, opens: string | undefined) {
  const color = ribbonsColor(o.hues === true);
  if (o.bands !== undefined) {
    s.track('bands', { label: 'Cytobands', height: 12 }).layer('bands', { type: 'bands', data: o.bands } as never);
    s.track('hues', { label: 'Chromosomes', height: 4 }).layer('hues', {
      type: 'bands',
      data: o.bands,
      color,
    } as never);
  }
  s.track('links', { label: 'Links', height: 4 }).layer('ribbons', {
    type: 'rondo:ribbons',
    fallback: 'pairs',
    data: o.links,
    color,
    ...(opens ? { opens } : {}),
  } as never);
}

/** Writes the genome ring into a view (_spec/README.md, The builder helper). */
export function rondoRing(view: ViewBuilder, o: RondoRingOptions): RondoRing {
  const id = o.id ?? 'ring';
  const detailId = `${id}-detail`;
  const below = o.detail === 'below';
  const opens = o.opens ? `${view.address}/${below ? detailId : id}` : undefined;
  const params = {
    of: o.assembly,
    layout: 'rings' as const,
    center: 'none' as const,
    ...(o.height !== undefined ? { height: o.height } : {}),
  };
  const ring = view.section(id, params as never, (s) => ringTracks(s, o, opens));
  const detail = below ? view.section(detailId, params as never, (s) => ringTracks(s, o, opens)) : undefined;
  if (o.square !== false) {
    const of = (detail ?? ring).address as Address<'section'>;
    view.navigator('square', {
      of: [of, of],
      axes: 'both',
      system: o.assembly,
      // its cells by count (`rondo:cells`); a reader without the extension draws the lattice (`axis`), said
      content: 'rondo:square',
      switch: 'include',
      size: [220, 220],
    } as never);
  }
  return { ring, ...(detail ? { detail } : {}) };
}
