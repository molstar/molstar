/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, isColumn, isNode, type Node, type Value } from '@molstar/core/repr/tree';
import { isColorSpelling } from '@molstar/seqx-core/theme/spelling';
import { ChannelWords, ClauseKinds, PaletteKinds, ThemeChannels, ThemeForms } from '@molstar/seqx-core/theme/forms';

// The checks of a layer's channel slots the description cannot make (_spec/README.md "Themes", Forms; the
// `seq.lang.theme.*` codes): a value of the channel's kind, colors that read, a scale that exists, edges and stops in
// order, comparisons within a type (against the column's declared type), nesting within `themeDepth`.

export type ThemeCode =
  | 'seq.lang.theme.form'
  | 'seq.lang.theme.color'
  | 'seq.lang.theme.channel'
  | 'seq.lang.theme.type'
  | 'seq.lang.theme.stops'
  | 'seq.lang.theme.scale'
  | 'seq.lang.theme.depth';

export interface ThemeCheckContext {
  /** Whether an address names a `scale` node. */
  isScale(address: string): boolean;
  /** A field's declared type on the layer's data: `number`, `text`, or undefined (unknown). */
  fieldType(name: string): 'number' | 'text' | undefined;
  readonly depth: number;
}

/** The declared type of each column of inline rows: a typed array holds numbers, a category or a list text. */
export function rowTypes(rows: Value | undefined): (name: string) => 'number' | 'text' | undefined {
  if (!rows || typeof rows !== 'object' || Array.isArray(rows) || isColumn(rows)) return () => undefined;
  const r = rows as Record<string, Value>;
  return (name) => {
    const c = r[name];
    if (c === undefined) return undefined;
    if (c instanceof CategoryColumn || Array.isArray(c)) return 'text';
    return isColumn(c) ? 'number' : undefined;
  };
}

/** The checks of one layer's channel slots: each issue with its code, message and the place in the layer's params. */
export function checkLayerThemes(
  layer: Node,
  ctx: ThemeCheckContext,
  report: (code: ThemeCode, message: string, ...tail: (string | number)[]) => void,
): void {
  for (const name of Object.keys(ThemeChannels) as (keyof typeof ThemeChannels)[]) {
    const v = layer.params?.[name];
    if (v === undefined) continue;
    const kind = ThemeChannels[name];
    if (name === 'ink' && v === 'auto') continue;
    const words = kind === 'word' ? ChannelWords[name as keyof typeof ChannelWords] : undefined;
    let deep = false;
    const value = (x: Value, depth: number, tail: (string | number)[]) => {
      if (isNode(x)) return form(x, depth, tail);
      if (typeof x === 'string') {
        if (kind === 'color' && !isColorSpelling(x))
          report('seq.lang.theme.color', `${JSON.stringify(x)} is no color`, ...tail);
        else if (kind === 'number') report('seq.lang.theme.channel', `The ${name} channel takes a number`, ...tail);
        else if (words && !words.includes(x))
          report('seq.lang.theme.channel', `The ${name} channel takes ${words.join(', ')}`, ...tail);
      } else if (typeof x === 'number' && kind !== 'number')
        report(
          'seq.lang.theme.channel',
          `The ${name} channel takes ${kind === 'color' ? 'a color' : 'a word'}`,
          ...tail,
        );
    };
    const form = (n: Node, depth: number, tail: (string | number)[]) => {
      if (depth > ctx.depth) {
        if (!deep) report('seq.lang.theme.depth', `Theme forms nest past ${ctx.depth} levels`, ...tail);
        deep = true;
        return;
      }
      const p = (n.params ?? {}) as Record<string, Value>;
      if (!(ThemeForms as readonly string[]).includes(n.kind)) {
        report('seq.lang.theme.form', `A ${n.kind} is no form of a channel`, ...tail);
        return;
      }
      switch (n.kind) {
        case 'uniform':
          if (p.value !== undefined) value(p.value, depth + 1, [...tail, 'value']);
          return;
        case 'field':
          field(p, tail);
          return;
        case 'rules':
          (n.args ?? []).forEach((r, i) => {
            if (!isNode(r)) return;
            const rp = (r.params ?? {}) as Record<string, Value>;
            const when = Array.isArray(rp.when) ? rp.when : [];
            when.forEach((c, k) => clause(c, [...tail, 'args', i, 'when', k]));
            if (rp.value !== undefined) value(rp.value, depth + 1, [...tail, 'args', i, 'value']);
          });
          if (p.else !== undefined) value(p.else, depth + 1, [...tail, 'else']);
          return;
        case 'letters':
          if (kind === 'color') {
            for (const [l, c] of Object.entries((p.colors ?? {}) as Record<string, Value>))
              if (typeof c === 'string' && !isColorSpelling(c))
                report('seq.lang.theme.color', `${JSON.stringify(c)} is no color`, ...tail, 'colors', l);
            (Array.isArray(p.groups) ? p.groups : []).forEach((g, i) => {
              const gv = isNode(g) ? g.params?.value : undefined;
              if (typeof gv === 'string' && !isColorSpelling(gv))
                report('seq.lang.theme.color', `${JSON.stringify(gv)} is no color`, ...tail, 'groups', i);
            });
          }
          return;
      }
    };
    const field = (p: Record<string, Value>, tail: (string | number)[]) => {
      if (typeof p.scale === 'string' && !ctx.isScale(p.scale))
        report('seq.lang.theme.scale', `${p.scale} names no scale`, ...tail, 'scale');
      const d = p.domain;
      if (Array.isArray(d) && typeof d[0] === 'number' && typeof d[1] === 'number' && d[0] > d[1])
        report('seq.lang.theme.stops', 'The domain’s ends are out of order', ...tail, 'domain');
      const pal = p.palette;
      if (pal === undefined) return;
      if (!isNode(pal) || !(PaletteKinds as readonly string[]).includes(pal.kind)) {
        report('seq.lang.theme.form', 'A palette goes here', ...tail, 'palette');
        return;
      }
      if (kind !== 'color')
        report('seq.lang.theme.channel', `The ${name} channel takes no palette (a scale’s range)`, ...tail, 'palette');
      palette(pal, [...tail, 'palette']);
    };
    const palette = (pal: Node, tail: (string | number)[]) => {
      const p = (pal.params ?? {}) as Record<string, Value>;
      const edges = p.edges;
      if (Array.isArray(edges) && !ascending(edges as number[]))
        report('seq.lang.theme.stops', 'The edges are out of order', ...tail, 'edges');
      const colors = p.colors;
      if (Array.isArray(colors)) {
        const stops = colors.every((c) => Array.isArray(c) && c.length === 2);
        if (stops && !ascending(colors.map((c) => (c as Value[])[1] as number)))
          report('seq.lang.theme.stops', 'The stops are out of order', ...tail, 'colors');
        colors.forEach((c, i) => {
          const s = Array.isArray(c) ? c[0] : c;
          if (typeof s === 'string' && !isColorSpelling(s))
            report('seq.lang.theme.color', `${JSON.stringify(s)} is no color`, ...tail, 'colors', i);
        });
      } else if (colors && typeof colors === 'object' && !isColumn(colors) && !isNode(colors))
        for (const [k, c] of Object.entries(colors as Record<string, Value>))
          if (typeof c === 'string' && !isColorSpelling(c))
            report('seq.lang.theme.color', `${JSON.stringify(c)} is no color`, ...tail, 'colors', k);
    };
    const clause = (c: Value, tail: (string | number)[]) => {
      if (!isNode(c) || !(ClauseKinds as readonly string[]).includes(c.kind)) {
        report('seq.lang.theme.form', 'A clause goes here', ...tail);
        return;
      }
      if (c.kind === 'in' || c.kind === 'within') return;
      const [m, ...rest] = c.args ?? [];
      const declared =
        isNode(m) && m.kind === 'field'
          ? ctx.fieldType(String(m.params?.name ?? ''))
          : isNode(m) && m.kind === 'value'
            ? 'number'
            : isNode(m) && m.kind === 'letter'
              ? 'text'
              : undefined;
      const numeric = c.kind === 'lt' || c.kind === 'le' || c.kind === 'gt' || c.kind === 'ge';
      for (const x of rest) {
        const t = typeof x === 'number' ? 'number' : typeof x === 'string' ? 'text' : undefined;
        if ((numeric && t !== 'number') || (declared && t && t !== declared))
          report(
            'seq.lang.theme.type',
            `The clause compares ${declared ?? 'a number'} with ${t ?? 'a value'}`,
            ...tail,
          );
      }
    };
    value(v, 0, [name]);
  }
}

function ascending(xs: readonly (number | null)[]): boolean {
  for (let i = 1; i < xs.length; i++) {
    const a = xs[i - 1];
    const b = xs[i];
    if (typeof a === 'number' && typeof b === 'number' && b < a) return false;
  }
  return true;
}
