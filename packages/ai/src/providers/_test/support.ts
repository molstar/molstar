/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type {
  JsonObject,
  Message,
  Provider,
  ModelFrame,
  ProjectedRequest,
  ProjectInput,
  ToolDefinition,
  ToolResultPart,
} from '../../messages/index.js';

// The fixtures replay history without provider signatures on purpose; the SDK warns about it.
(globalThis as { AI_SDK_LOG_WARNINGS?: boolean }).AI_SDK_LOG_WARNINGS = false;

export const lookupTool: ToolDefinition = {
  tool: { id: 'fixture.lookup' },
  name: 'fixture.lookup',
  description: 'Look up one fixture value.',
  inputSchema: {
    type: 'object',
    properties: { id: { type: 'string' } },
    required: ['id'],
  },
};

export type Vendor = 'anthropic' | 'google' | 'openai';
export type InputOptions = Partial<Omit<ProjectInput, 'capabilities'>>;

/** A project input without its capabilities, which `project` resolves from the provider. */
export function projectInput(
  provider: Vendor,
  model: string,
  extras: InputOptions = {},
): Omit<ProjectInput, 'capabilities'> {
  return {
    configuration: { target: { provider, model }, parameters: { maxOutputTokens: 128 } },
    messages: [{ role: 'user', parts: [{ kind: 'text', text: 'Look up 1abc.' }] }],
    tools: [lookupTool],
    ...extras,
  };
}

/** Project with the capabilities the provider reports for the input's configuration, as the harness does. */
export async function project(
  provider: Provider,
  input: Omit<ProjectInput, 'capabilities'>,
  capabilities?: ProjectInput['capabilities'],
): Promise<ProjectedRequest> {
  return await provider.project({
    ...input,
    capabilities: capabilities ?? (await provider.capabilities(input.configuration)),
  });
}

export const payloadOf = (request: ProjectedRequest): JsonObject => request.payload as JsonObject;

export function resultOf(frames: readonly ModelFrame[]): Extract<ModelFrame, { kind: 'result' }> {
  const result = frames.find((frame) => frame.kind === 'result');
  if (result?.kind !== 'result') throw new Error('Expected a result frame');
  return result;
}

/** A tool result answering the call. */
export const resultPart = (requestId: string, text = 'found'): ToolResultPart => ({
  kind: 'tool-result',
  requestId,
  tool: lookupTool.tool,
  status: 'succeeded',
  parts: [{ kind: 'text', text }],
});

/** A user message answering the call. */
export const answer = (requestId: string, text = 'found'): Message => ({
  role: 'user',
  parts: [resultPart(requestId, text)],
});

/** A system message declaring a tool: `deferred: false` also loads it, where the message stands. */
export const declares = (definition: ToolDefinition, deferred = false): Message => ({
  role: 'system',
  parts: [{ kind: 'tool-definition', definition, deferred }],
});

export async function collectFrames(
  provider: Provider,
  request: ProjectedRequest,
  signal: AbortSignal = new AbortController().signal,
): Promise<readonly ModelFrame[]> {
  const frames: ModelFrame[] = [];
  for await (const frame of provider.invoke({ request, signal })) {
    frames.push(frame);
  }
  return frames;
}

export function anthropicSse(text: string, tool?: { readonly name: string; readonly input: string }): string {
  const events: readonly (readonly [string, unknown])[] = [
    [
      'message_start',
      {
        type: 'message_start',
        message: {
          id: 'msg_fixture',
          type: 'message',
          role: 'assistant',
          content: [],
          model: 'claude-sonnet-4-5',
          stop_reason: null,
          usage: { input_tokens: 12, output_tokens: 0, cache_read_input_tokens: 4, cache_creation_input_tokens: 1 },
        },
      },
    ],
    ['content_block_start', { type: 'content_block_start', index: 0, content_block: { type: 'text', text: '' } }],
    ['content_block_delta', { type: 'content_block_delta', index: 0, delta: { type: 'text_delta', text } }],
    ['content_block_stop', { type: 'content_block_stop', index: 0 }],
    ...(tool === undefined
      ? []
      : ([
          [
            'content_block_start',
            {
              type: 'content_block_start',
              index: 1,
              content_block: { type: 'tool_use', id: 'toolu_fixture', name: tool.name, input: {} },
            },
          ],
          [
            'content_block_delta',
            { type: 'content_block_delta', index: 1, delta: { type: 'input_json_delta', partial_json: tool.input } },
          ],
          ['content_block_stop', { type: 'content_block_stop', index: 1 }],
        ] as const)),
    [
      'message_delta',
      {
        type: 'message_delta',
        delta: { stop_reason: tool === undefined ? 'end_turn' : 'tool_use' },
        usage: { output_tokens: 11, cache_read_input_tokens: 4, cache_creation_input_tokens: 1 },
      },
    ],
    ['message_stop', { type: 'message_stop' }],
  ];
  return events.map(([event, data]) => `event: ${event}\ndata: ${JSON.stringify(data)}\n\n`).join('');
}

export function anthropicThinkingSse(
  thinking: { readonly text: string; readonly signature: string } | { readonly redactedData: string },
  tool: { readonly name: string; readonly input: string },
): string {
  const events: readonly (readonly [string, unknown])[] = [
    [
      'message_start',
      {
        type: 'message_start',
        message: {
          id: 'msg_fixture',
          type: 'message',
          role: 'assistant',
          content: [],
          model: 'claude-sonnet-5',
          stop_reason: null,
          usage: { input_tokens: 12, output_tokens: 0 },
        },
      },
    ],
    ...('redactedData' in thinking
      ? ([
          [
            'content_block_start',
            {
              type: 'content_block_start',
              index: 0,
              content_block: { type: 'redacted_thinking', data: thinking.redactedData },
            },
          ],
        ] as const)
      : ([
          [
            'content_block_start',
            { type: 'content_block_start', index: 0, content_block: { type: 'thinking', thinking: '' } },
          ],
          [
            'content_block_delta',
            { type: 'content_block_delta', index: 0, delta: { type: 'thinking_delta', thinking: thinking.text } },
          ],
          [
            'content_block_delta',
            {
              type: 'content_block_delta',
              index: 0,
              delta: { type: 'signature_delta', signature: thinking.signature },
            },
          ],
        ] as const)),
    ['content_block_stop', { type: 'content_block_stop', index: 0 }],
    [
      'content_block_start',
      {
        type: 'content_block_start',
        index: 1,
        content_block: { type: 'tool_use', id: 'toolu_fixture', name: tool.name, input: {} },
      },
    ],
    [
      'content_block_delta',
      { type: 'content_block_delta', index: 1, delta: { type: 'input_json_delta', partial_json: tool.input } },
    ],
    ['content_block_stop', { type: 'content_block_stop', index: 1 }],
    ['message_delta', { type: 'message_delta', delta: { stop_reason: 'tool_use' }, usage: { output_tokens: 11 } }],
    ['message_stop', { type: 'message_stop' }],
  ];
  return events.map(([event, data]) => `event: ${event}\ndata: ${JSON.stringify(data)}\n\n`).join('');
}

export function googleSse(
  text: string,
  tool?: { readonly name: string; readonly args: Record<string, string>; readonly thoughtSignature?: string },
): string {
  return `data: ${JSON.stringify({
    candidates: [
      {
        content: {
          role: 'model',
          parts: [
            { text },
            ...(tool === undefined
              ? []
              : [
                  {
                    functionCall: { id: 'call_fixture', name: tool.name, args: tool.args },
                    ...(tool.thoughtSignature === undefined ? {} : { thoughtSignature: tool.thoughtSignature }),
                  },
                ]),
          ],
        },
        finishReason: 'STOP',
      },
    ],
    usageMetadata: {
      promptTokenCount: 12,
      candidatesTokenCount: 11,
      cachedContentTokenCount: 4,
      thoughtsTokenCount: 3,
    },
  })}\n\n`;
}

export function recordingFetch(body: string): {
  readonly fetch: typeof globalThis.fetch;
  readonly requests: Array<{ readonly url: string; readonly headers: Headers; readonly body: string }>;
} {
  const requests: Array<{ readonly url: string; readonly headers: Headers; readonly body: string }> = [];
  return {
    requests,
    fetch: async (input, init) => {
      requests.push({
        url: String(input),
        headers: new Headers(init?.headers),
        body: typeof init?.body === 'string' ? init.body : '',
      });
      return new Response(body, { status: 200, headers: { 'content-type': 'text/event-stream' } });
    },
  };
}

export function anthropicErrorSse(): string {
  return `event: error\ndata: ${JSON.stringify({
    type: 'error',
    error: { type: 'api_error', message: 'provider failed' },
  })}\n\n`;
}

export function hangingFetch(): {
  readonly fetch: typeof globalThis.fetch;
  readonly started: Promise<void>;
} {
  let markStarted: (() => void) | undefined;
  const started = new Promise<void>((resolve) => {
    markStarted = resolve;
  });
  return {
    started,
    fetch: async (_input, init) => {
      markStarted?.();
      const signal = init?.signal;
      return new Promise((_, reject) => {
        const abort = () => {
          reject(signal?.reason instanceof Error ? signal.reason : new Error('aborted'));
        };
        if (signal?.aborted === true) {
          abort();
          return;
        }
        signal?.addEventListener('abort', abort, { once: true });
      });
    },
  };
}

export function openaiSse(
  text: string,
  tool?: { readonly name: string; readonly input: string },
  reasoning?: { readonly text: string; readonly encryptedContent: string },
): string {
  const events: unknown[] = [
    { type: 'response.created', response: { id: 'resp_fixture', created_at: 1, model: 'gpt-6-luna' } },
  ];
  let index = 0;
  if (reasoning !== undefined) {
    events.push(
      {
        type: 'response.output_item.added',
        output_index: index,
        item: { type: 'reasoning', id: 'rs_fixture', encrypted_content: null },
      },
      {
        type: 'response.reasoning_summary_part.added',
        output_index: index,
        item_id: 'rs_fixture',
        summary_index: 0,
        part: { type: 'summary_text', text: '' },
      },
      {
        type: 'response.reasoning_summary_text.delta',
        output_index: index,
        item_id: 'rs_fixture',
        summary_index: 0,
        delta: reasoning.text,
      },
      {
        type: 'response.output_item.done',
        output_index: index,
        item: { type: 'reasoning', id: 'rs_fixture', encrypted_content: reasoning.encryptedContent },
      },
    );
    index++;
  }
  if (text.length > 0) {
    events.push(
      { type: 'response.output_item.added', output_index: index, item: { type: 'message', id: 'msg_fixture' } },
      { type: 'response.output_text.delta', output_index: index, item_id: 'msg_fixture', delta: text },
      { type: 'response.output_item.done', output_index: index, item: { type: 'message', id: 'msg_fixture' } },
    );
    index++;
  }
  if (tool !== undefined) {
    events.push(
      {
        type: 'response.output_item.added',
        output_index: index,
        item: { type: 'function_call', id: 'fc_fixture', call_id: 'call_fixture', name: tool.name, arguments: '' },
      },
      { type: 'response.function_call_arguments.delta', output_index: index, item_id: 'fc_fixture', delta: tool.input },
      {
        type: 'response.output_item.done',
        output_index: index,
        item: {
          type: 'function_call',
          id: 'fc_fixture',
          call_id: 'call_fixture',
          name: tool.name,
          arguments: tool.input,
          status: 'completed',
        },
      },
    );
  }
  events.push({
    type: 'response.completed',
    response: {
      usage: {
        input_tokens: 12,
        input_tokens_details: { cached_tokens: 4 },
        output_tokens: 11,
        output_tokens_details: { reasoning_tokens: 3 },
      },
    },
  });
  return events.map((event) => `data: ${JSON.stringify(event)}\n\n`).join('');
}
