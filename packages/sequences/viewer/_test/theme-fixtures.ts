/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Theme fixtures (LY3a): a theme context over a hand-made layer through the default registry, a reused location.

import type { ChannelEntry, ChannelKind, ResolvedLayer, ResolvedSection, ThemeRef } from '@molstar/seq-core/format';
import type { Issue, KeyTarget, SeqSet } from '@molstar/seq-core/model';
import { createSeqViewer, type Theme, ThemeLocation } from '../src/index.js';

import { createTheme, themeContext } from '../src/themes.js';

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const idOf = (set: SeqSet, key: string) => (set.byName(key, t) ? t.id : -1);
export const channel = (kind: ChannelKind, unit: ChannelEntry['unit'] = null): ChannelEntry => ({
  name: kind,
  kind,
  unit,
  default: null,
});

interface Opts {
  data?: string[];
  blocks?: { unit?: string; alignment?: string }[];
  dark?: boolean;
  domain?: [number, number];
  /** The layer's type (what a color scheme fits reads); none by default. */
  type?: string;
}
/** A theme context over a hand-made layer of `set`, through the default registry. */
export function context(set: SeqSet, o: Opts = {}) {
  const registry = createSeqViewer().registry;
  const issues: Issue[] = [];
  const data = (o.data ?? []).map((k) => idOf(set, k));
  const [domainMin, domainMax] = o.domain ?? [Number.NaN, Number.NaN];
  const layer = { ref: 'L', data, kept: data.map(() => null), domainMin, domainMax, params: {}, placements: [] };
  if (o.type) Object.assign(layer, { type: set.strings.intern(o.type) });
  const blocks = (o.blocks ?? [{ unit: 'A' }]).map((b) => ({
    unit: b.unit ? idOf(set, b.unit) : -1,
    alignment: b.alignment ? idOf(set, b.alignment) : -1,
  }));
  const input = {
    set,
    section: { blocks } as unknown as ResolvedSection,
    layer: layer as unknown as ResolvedLayer,
    annotations: data.map((id) => set.annotations[id] as NonNullable<SeqSet['annotations'][number]>),
  };
  const report = (i: Issue) => issues.push(i);
  const tc = themeContext(input, o.dark ?? false, registry, report);
  const make =
    (kind: ChannelKind) =>
    (theme: string, params: Record<string, unknown> = {}, fallback?: ThemeRef) =>
      createTheme(
        tc,
        { theme, params, fallback: fallback ?? null },
        channel(kind, kind === 'size' ? 'px' : null),
        registry,
        report,
      ) as Theme;
  return { tc, issues, registry, color: make('color'), glyph: make('glyph'), size: make('size') };
}
const loc = ThemeLocation.create();
export const at = (o: Partial<ThemeLocation>) => {
  ThemeLocation.reset(loc);
  return Object.assign(loc, o);
};
