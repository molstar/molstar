/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { JsonObject, ProjectedRequest } from '../../messages/index.js';
import { decodeProjectedPayload } from '../common/payload.js';

/** Private SDK conversion placeholder, bound to one retained request and item index.
 * The installed SDK has no additional_tools input type. The transport replaces
 * validated developer placeholders before HTTP, leaving all other SDK items intact. */
export function openAIDisclosureMarker(requestDigest: string, index: number): string {
  return `molstar-ai:openai:additional-tools:${requestDigest}:${index}`;
}

export function openAIDisclosureTransport(
  request: ProjectedRequest,
  fetch: typeof globalThis.fetch,
): typeof globalThis.fetch {
  const payload = decodeProjectedPayload(request.payload);
  if (payload.appendsTools !== true) return fetch;
  const additions = payload.messages.flatMap((message, index) =>
    message.type === 'additional_tools' && Array.isArray(message.tools)
      ? [{ marker: openAIDisclosureMarker(request.digest, index), tools: message.tools }]
      : [],
  );
  return async (url, init) => {
    if (typeof init?.body !== 'string') throw new Error('OpenAI disclosure transport expected an SDK JSON request.');
    const body = JSON.parse(init.body) as Record<string, unknown>;
    if (!Array.isArray(body.input)) throw new Error('OpenAI disclosure transport expected ordered input items.');
    let next = 0;
    body.input = body.input.map((entry) => {
      if (entry === null || typeof entry !== 'object' || Array.isArray(entry)) return entry;
      const item = entry as JsonObject;
      const candidate = additions.find((addition) => addition.marker === item.content);
      if (candidate === undefined) return entry;
      if (item.role !== 'developer' || candidate !== additions[next])
        throw new Error('OpenAI disclosure placeholder order or role changed during SDK conversion.');
      next++;
      return { type: 'additional_tools', role: 'developer', tools: candidate.tools };
    });
    if (next !== additions.length) throw new Error('OpenAI disclosure placeholder was lost during SDK conversion.');
    // Tools were supplied to the SDK solely to parse returned calls. Native full
    // definitions already occupy their exact positions in the input above.
    delete body.tools;
    return fetch(url, { ...init, body: JSON.stringify(body) });
  };
}
