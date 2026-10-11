/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A text-token profile: ASCII characters per token, and UTF-8 bytes per token for other text. */
export interface TokenEstimate {
  readonly asciiCharactersPerToken: number;
  readonly otherBytesPerToken: number;
}

/** The uncalibrated profile: ASCII prose is about four characters a token. */
export const DEFAULT_ESTIMATE: TokenEstimate = { asciiCharactersPerToken: 4, otherBytesPerToken: 2 };

/** An offline text-token heuristic, not a tokenizer; an application calibrates it per provider
 * (`tokenEstimate` on the adapter's options). Only `inputTokens` uses it. */
export function estimateTextTokens(text: string, profile: TokenEstimate = DEFAULT_ESTIMATE): number {
  let ascii = 0;
  let other = 0;
  for (let i = 0; i < text.length; i += 1) {
    const code = text.charCodeAt(i);
    if (code < 0x80) ascii += 1;
    else if (code < 0x800) other += 2;
    // A surrogate pair is one 4-byte character: 2 for each half.
    else other += code >= 0xd800 && code <= 0xdfff ? 2 : 3;
  }
  return Math.ceil(ascii / profile.asciiCharactersPerToken + other / profile.otherBytesPerToken);
}
