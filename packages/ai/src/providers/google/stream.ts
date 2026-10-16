/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StreamPolicy } from '../common/stream.js';

export const streamPolicy: StreamPolicy = {
  options(payload) {
    return {
      ...(payload.system.length > 0 ? { system: payload.system } : {}),
      ...(payload.reasoning !== undefined
        ? {
            providerOptions: {
              google: {
                thinkingConfig:
                  payload.reasoning.mode === 'budget'
                    ? { thinkingBudget: payload.reasoning.tokens }
                    : { ...(payload.reasoning.effort ? { thinkingLevel: payload.reasoning.effort } : {}) },
              },
            },
          }
        : {}),
    };
  },
  inlineMessage(message, _index, _digest) {
    if (message.role === 'system') throw new Error('Projected provider does not support inline system messages.');
    return undefined;
  },
  replayOptions(value) {
    const options = googleReplayOptions(value);
    return options === undefined ? undefined : { google: { ...options } };
  },
  createReplay() {
    return {
      reasoning() {},
      toolCall(call, part) {
        const thoughtSignature = part.providerMetadata?.google?.thoughtSignature;
        return typeof thoughtSignature === 'string' && thoughtSignature.length > 0
          ? { ...call, providerOptions: { google: { thoughtSignature } } }
          : call;
      },
      finish: (calls) => calls,
    };
  },
};

function googleReplayOptions(value: unknown): { readonly thoughtSignature: string } | undefined {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) {
    return undefined;
  }
  const google = (value as Record<string, unknown>).google;
  if (google === null || typeof google !== 'object' || Array.isArray(google)) {
    return undefined;
  }
  const thoughtSignature = (google as Record<string, unknown>).thoughtSignature;
  return typeof thoughtSignature === 'string' && thoughtSignature.length > 0 ? { thoughtSignature } : undefined;
}
