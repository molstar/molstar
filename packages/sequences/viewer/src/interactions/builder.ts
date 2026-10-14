/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Builder functions (seqspec:libraries#interactions): the `computed` fingerprint node, and theme references whose
// fallbacks write this vocabulary's colors and glyphs out as kernel themes, so a reader without the extension draws the
// same.

import type { DataBuilder, NodeHandle } from "@molstar/seq-core/builder";
import type { ThemeRef } from "@molstar/seq-core/format";
import { type KindName, KindStyle, Kinds } from "./vocabulary.js";

export type InteractionRef = string | NodeHandle;
export interface FingerprintParams {
  ref?: string;
  /** Pair annotations. */
  pairs: InteractionRef[];
  /** Units: ligands or poses, a row each. */
  ligands: InteractionRef[];
  /** A unit, a sequence or an alignment: the columns' coordinate system. */
  system: InteractionRef;
  /** Only these kinds; null: all. */
  kinds?: KindName[] | null;
  custom?: Record<string, unknown>;
}

/** The `interactions:fingerprint` computed node: a matrix over the system, one state per ligand. */
export function fingerprint(data: DataBuilder, params: FingerprintParams): NodeHandle {
  const { ref, custom, ...options } = params;
  return data.computed({ ref, custom, provider: "interactions:fingerprint", options });
}

const per = <T>(f: (k: KindName) => T) => Object.fromEntries(Kinds.map((k) => [k, f(k)]));

/** The `interactions:kind` color theme; its fallback a `field_palette` over `category`. */
export function interactionColors(): ThemeRef {
  const colors = per((k) => KindStyle[k][1]);
  const palette = { kind: "categorical", colors, missing_color: KindStyle.other[1] };
  return {
    theme: "interactions:kind",
    params: {},
    fallback: { theme: "field_palette", params: { field: "category", palette } },
  };
}

/** The `interactions:kind` glyph theme; its fallback a `field_map` over `category`. */
export function interactionGlyphs(): ThemeRef {
  const glyph = (k: KindName) => ({ shape: KindStyle[k][2], label: KindStyle[k][0] });
  const params = { field: "category", map: per(glyph), default: glyph("other") };
  return { theme: "interactions:kind", params: {}, fallback: { theme: "field_map", params } };
}
