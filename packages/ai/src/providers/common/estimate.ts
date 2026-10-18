/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { base64Bytes, estimateAttachmentTokens } from './attachments.js';

/** A text-token profile: ASCII characters per token, and UTF-8 bytes per token for other text. */
export interface TokenEstimate {
  readonly asciiCharactersPerToken: number;
  readonly otherBytesPerToken: number;
}

/** The uncalibrated profile: ASCII prose is about four characters a token. */
export const DEFAULT_ESTIMATE: TokenEstimate = { asciiCharactersPerToken: 4, otherBytesPerToken: 2 };

/** An offline text-token heuristic, not a tokenizer; an application calibrates it per provider
 * (`tokenEstimate` on the provider's options). Only `inputTokens` uses it. */
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

/** A tool's schema text plus a fixed allowance. */
export function estimateTool(
  tool: { readonly name: string; readonly description: string; readonly inputSchema: unknown },
  profile?: TokenEstimate,
): number {
  return 16 + estimateTextTokens(tool.name + tool.description + JSON.stringify(tool.inputSchema), profile);
}

type Part = { readonly [key: string]: unknown };

/** One part of a message's content: text and replayed thinking count, files by their modality, signatures not. */
export function estimatePart(part: Part, text: (value: string) => number): number {
  switch (part.type) {
    case 'text':
    case 'reasoning':
      return 4 + text(String(part.text ?? ''));
    case 'file': {
      const data = part.data as Part;
      return (
        4 +
        (data.type === 'text'
          ? text(String(data.text ?? ''))
          : estimateAttachmentTokens(String(part.mediaType ?? ''), base64Bytes(String(data.data ?? ''))))
      );
    }
    case 'tool-call':
      return 8 + text(JSON.stringify(part.input ?? {})) + text(String(part.toolName ?? ''));
    case 'tool-result': {
      const output = part.output as Part;
      return (
        8 +
        (Array.isArray(output.value)
          ? output.value.reduce<number>((sum, item) => sum + estimatePart(item as Part, text), 0)
          : text(typeof output.value === 'string' ? output.value : ''))
      );
    }
    default:
      return 0;
  }
}

/** A message's content (a string, or parts). */
export function estimateContent(content: unknown, text: (value: string) => number): number {
  if (typeof content === 'string') return text(content);
  return Array.isArray(content) ? content.reduce<number>((sum, part) => sum + estimatePart(part as Part, text), 0) : 0;
}
