/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The theme sublanguage's forms as data (_spec/README.md "Themes", Forms): a slot's body is a string (the `uniform`
// form), a number (a uniform size or opacity), or a node of `seq.theme` without an id.

/** A node of the theme sublanguage, as the tree holds it. */
export interface ThemeNode {
  readonly kind: string;
  readonly params?: Readonly<Record<string, unknown>>;
  readonly args?: readonly unknown[];
}

/** A channel slot's body. */
export type ThemeBody = string | number | ThemeNode;

/** The channels a layer may read, and what each holds. */
export const ThemeChannels = {
  color: 'color',
  ink: 'color',
  size: 'number',
  opacity: 'number',
  line: 'word',
  arrow: 'word',
  shape: 'word',
} as const;

export type ChannelName = keyof typeof ThemeChannels;
export type ChannelKind = (typeof ThemeChannels)[ChannelName];

/** The words a word channel takes. */
export const ChannelWords: Readonly<Record<'line' | 'arrow' | 'shape', readonly string[]>> = {
  line: ['solid', 'dashed', 'dotted', 'wavy', 'double', 'box', 'highlight'],
  arrow: ['none', 'end', 'start', 'both'],
  shape: ['circle', 'square', 'triangle', 'diamond'],
};

/** The forms (kinds a slot's node may be). */
export const ThemeForms = ['uniform', 'field', 'rules', 'scheme', 'letters'] as const;
export type ThemeFormKind = (typeof ThemeForms)[number];

/** The palettes (a `field` form's or a scale's `palette`). */
export const PaletteKinds = ['categorical', 'direct', 'binned', 'continuous'] as const;
export type PaletteKind = (typeof PaletteKinds)[number];

/** The clauses of a rule (data forms: a measure, then values); the words spell `=`, `!=`, `<`, `<=`, `>`, `>=`. */
export const ClauseKinds = ['eq', 'ne', 'lt', 'le', 'gt', 'ge', 'one-of', 'in', 'within'] as const;
export type ClauseKind = (typeof ClauseKinds)[number];

/** The built-in measures a clause reads (`field` takes `:name`); the column measures are registered (`share`, ...). */
export const BuiltInMeasures = ['field', 'value', 'letter'] as const;

/** A scale's params, as a `scale` node or a `field` form holds them inline. */
export interface ScaleParams {
  readonly type?: 'linear' | 'log' | 'sqrt';
  readonly domain?: readonly [number | null, number | null];
  readonly palette?: ThemeNode;
  readonly range?: readonly [number, number];
  readonly mid?: number;
  readonly nice?: boolean;
  readonly label?: string;
  readonly unit?: string;
  readonly non_positive?: 'missing' | 'clamp';
}

export function isThemeNode(v: unknown): v is ThemeNode {
  return typeof v === 'object' && v !== null && !Array.isArray(v) && typeof (v as ThemeNode).kind === 'string';
}

/** A param of a node, or undefined. */
export function param(node: ThemeNode, name: string): unknown {
  const p = node.params;
  return p && Object.hasOwn(p, name) ? p[name] : undefined;
}
