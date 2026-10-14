/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `DrawResources` over a viewer context (seqspec:viewer#the-render-backend): a FontId's CSS font and size from the text
// manager, a skin token's color from the context's skin, a glyph from the registry's glyph table. The SVG backend
// resolves ids through it; the Canvas2D surface wraps it with a font size cache by id (its hot path).

import type { SeqViewerContext } from "../context.ts";
import type { Skin } from "../graph.ts";
import { fontPx } from "../text.ts";
import { glyphTableOf } from "../themes.ts";
import { type DrawResources, SkinToken } from "./backend.ts";

/** Skin tokens by `-token - 2` (`SkinToken`), as `Skin` keys. */
const TOKENS: (keyof Skin)[] = [];
for (const [name, id] of Object.entries(SkinToken)) TOKENS[-id - 2] = name as keyof Skin;

/** The resources of `ctx`, reading its skin and fonts as they stand when asked (a skin change shows at once). */
export function drawResources(ctx: SeqViewerContext): DrawResources {
  const glyphs = glyphTableOf(ctx.registry);
  return {
    font: (id) => ctx.text.fontString(id),
    fontPx: (id) => fontPx(ctx.text.fontString(id)),
    token: (t) => ctx.skin.value[TOKENS[-t - 2]] as string,
    glyph: (id) => glyphs.glyph(id),
  };
}
