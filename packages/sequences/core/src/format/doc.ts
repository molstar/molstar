/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The document's contract types (seqspec:sources#the-loader): the envelope, nodes, snapshots, patches, invocations,
// theme references.

/** The format version this package reads and writes (seqspec:format#versioning). */
export const SvsVersion = '1.0.0';
export const SvsFormat = 'seq_view_spec';

export interface DocNode {
  kind: string;
  ref?: string;
  params?: Record<string, unknown>;
  custom?: Record<string, unknown>;
  children?: DocNode[];
}

export interface Metadata {
  title?: string | null;
  description?: string | null;
  description_format?: 'markdown' | 'plaintext';
  timestamp?: string | null;
  authors?: string[] | null;
}

export interface SeqViewSpecDoc {
  format: 'seq_view_spec';
  version: string;
  revision?: string | null;
  metadata?: Metadata;
  requires?: string[];
  root: DocNode;
  snapshots?: Snapshot[];
}

export interface Snapshot {
  ref: string;
  title?: string | null;
  description?: string | null;
  description_format?: 'markdown' | 'plaintext';
  ops?: Op[];
  duration_ms?: number;
  transition_ms?: number;
}

export interface Patch {
  format: 'seq_view_spec_patch';
  base?: string;
  revision?: string;
  ops: Op[];
}

/** A sibling's ref, `"first"` or `"last"`. */
export type Anchor = string;
export type Op =
  | { op: 'set'; ref: string; params?: Record<string, unknown>; custom?: Record<string, unknown> }
  | { op: 'remove'; ref: string }
  | { op: 'insert'; parent: string; after?: Anchor; node: DocNode }
  | { op: 'move'; ref: string; parent: string; after?: Anchor };

/**
 * A request to a viewer, plain JSON (seqspec:format#descriptions-and-invocations): a registered command's id and its
 * params (`{}` by default).
 */
export interface Invocation {
  readonly command: string;
  readonly params?: Readonly<Record<string, unknown>>;
}

export interface ThemeRef {
  theme: string;
  params?: Record<string, unknown>;
  fallback?: ThemeRef | null;
  /** Set by resolution on a named rule table's expansion (seqspec:format#themes): the table it stands for. */
  table?: string;
}

export type GlyphShape = 'circle' | 'diamond' | 'triangle' | 'square' | 'pin' | 'star' | 'cross' | 'tick';
/** How a range's item is drawn (seqspec:format#themes): `solid` (none given) is the plain bar or underline. */
export type LineStyle = 'solid' | 'dashed' | 'dotted' | 'wavy' | 'double' | 'box' | 'highlight';
export type Glyph =
  | GlyphShape
  | { shape?: GlyphShape; char?: string; path?: string; box?: [number, number]; label?: string; fallback?: GlyphShape };

/** The model's, not a second declaration of it. */
export type { Geometry } from '@molstar/seq-core/model';
