/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LegendItem, LegendRamp } from '@molstar/prism-ui';
import { type ColorSchemeEntry, type ColorStop, KernelColorSchemes, schemeLegend } from '@molstar/seq-core/format';

// What the colors of the letters mean (prism-sequence.md §6): read from the color scheme the model chose, its public
// definition (`KernelColorSchemes`), in Prism's `Legend` terms. Seq*'s own legend model stays inside its UIs.

export interface SequenceLegend {
  /** The scheme's name as a picker shows it; absent: the letters are not colored. */
  readonly title?: string;
  readonly items: readonly LegendItem[];
  readonly ramps: readonly LegendRamp[];
}

export const NoLegend: SequenceLegend = { items: [], ramps: [] };

/** The registered scheme `name` (the kernel's). */
export function colorScheme(name: string): ColorSchemeEntry | undefined {
  return KernelColorSchemes.find((e) => e.name === name);
}

const stopsOf = (colors: readonly ColorStop[]): { t: number; color: string }[] =>
  colors.map((c, i) => ({
    t: typeof c === 'string' ? (colors.length === 1 ? 0 : i / (colors.length - 1)) : c[1],
    color: typeof c === 'string' ? c : c[0],
  }));

/** The legend of a scheme for a light or dark ink; empty for no scheme or an unknown name. */
export function legendOf(scheme: string | null, dark: boolean): SequenceLegend {
  const entry = scheme === null ? undefined : colorScheme(scheme);
  if (!entry) return NoLegend;
  const def = entry.scheme;
  const title = entry.label;
  if (def.kind === 'letters' || def.kind === 'theme') {
    const items = schemeLegend(undefined, { custom: def }, dark).map((i) => ({ label: i.label, color: i.color }));
    return { title, items, ramps: [] };
  }
  if (def.kind === 'categorical') {
    const colors = dark && def.dark && def.dark !== 'same' ? def.dark : def.colors;
    const items = Array.isArray(colors)
      ? colors.map((color, i) => ({ label: `${i + 1}`, color }))
      : Object.entries(colors).map(([label, color]) => ({ label, color: color as string }));
    return { title, items, ramps: [] };
  }
  const colors = dark && def.dark && def.dark !== 'same' ? def.dark : def.colors;
  return { title, items: [], ramps: [{ stops: stopsOf(colors), lo: 'Low', hi: 'High' }] };
}
