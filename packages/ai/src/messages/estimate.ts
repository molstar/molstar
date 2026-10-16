/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * A text's size in tokens, plain: about four characters a token, at least one. The trace shows it; the providers
 * calibrate their own admission estimate (`providers/common/estimate.ts`), and compaction sizes with a cautious three
 * characters a token.
 */
export const estimatedTokens = (text: string): number => Math.max(1, Math.round(text.length / 4));
