/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json, ModelVisiblePart, ToolDefinition } from '@molstar/ai-agent';
import { HarnessError } from '../contracts/ports.js';
import type { ToolFrame, ToolHandler } from '../contracts/runner.js';

const outside = () => {
  throw new HarnessError('invalid', 'a tool called without a harness has no journal and no nested runs');
};

/**
 * A tool called without a model (harness.md §8): no harness, the caller's key, the same frame protocol.
 * Its journal reads nothing and refuses writes; `run()` is refused.
 */
export async function callTool(
  handler: ToolHandler,
  input: {
    readonly definition: ToolDefinition;
    readonly arguments: Json;
    readonly key: string;
    readonly signal?: AbortSignal;
  },
): Promise<{ status: 'succeeded' | 'failed'; parts: ModelVisiblePart[] }> {
  const failed = (text: string) => ({ status: 'failed' as const, parts: [{ kind: 'text' as const, text }] });
  let result: Extract<ToolFrame, { kind: 'result' }> | undefined;
  try {
    for await (const frame of handler.invoke({
      harnessId: '',
      turnId: '',
      callId: input.key,
      path: [],
      key: input.key,
      prior: 'none',
      attempt: 1,
      definition: input.definition,
      arguments: input.arguments,
      signal: input.signal ?? new AbortController().signal,
      journal: { read: async () => ({ messages: [] }), append: outside, edit: outside },
      run: async () => outside(),
    })) {
      if (result !== undefined) return failed('stream-protocol: a frame after the result');
      if (frame.kind === 'result') result = frame;
    }
  } catch (error) {
    return failed(String(error));
  }
  return result === undefined
    ? failed('stream-protocol: no result')
    : { status: result.status, parts: [...result.parts] };
}
