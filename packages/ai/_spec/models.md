# Messages and models (`src/messages/`, `src/providers/` of `@molstar/ai`)

## Messages (`src/messages/messages.ts`)

- `Message` is `SystemMessage | UserMessage | AssistantMessage`, each `{ role, parts, label? }`. Parts:
  - `ModelVisiblePart`: `{ kind: 'text', text }` or `{ kind: 'reference', reference: Reference }` (a URI with
    `mediaType`, `byteLength`, `digest`; providers inline only `data:` URIs, others become text `[ref <uri>]`).
  - system only: `ToolDefinitionPart { definition, deferred }`; `deferred: false` loads the tool, `true` only declares
    it. Never sent as text.
  - user only: `ToolResultPart { requestId, tool, status, parts }`.
  - assistant only: `ToolCallPart { requestId, tool, arguments, providerOptions? }` (`providerOptions` keeps the
    provider's replay data: thinking signatures, encrypted reasoning).
- `ToolDefinition { tool: { id }, name, description, inputSchema }`; the provider-safe wire name is derived from the id.
- `hashMessage(message)` (MurmurHash3_x86_128 over a structural encoding, synchronous) is an entry's id; `digest(json)`
  addresses blobs and heads; `canonicalJson`, `encodeMessage` are the encodings (`src/messages/hash.ts`).
- **`label?: EntryLabel { kind, key?, data? }`** says what an entry is, for the code that reads the journal: `kind` is
  namespaced by the package that writes it, `key` names which one, `data` holds small facts the parts also say.
  - Part of the message, so of its hash (an absent label leaves the hash as before; the same text under another label is
    another entry).
  - Never sent: providers read only `role` and `parts` (`src/providers/_test/labels.test.ts`), and sizes never count it.
    The harness never reads it.
  - Only what is true of the content; what happened when (turn, step, time, cost) is the event log's. There is no other
    metadata on a message.
  - The kinds in use: `session:request` (data `{ prompt, focus? }`), `session:guidance`, `session:guidance-update` (key:
    the text's fingerprint), `session:carry-over`, `session:memory` (`session/` `LABEL_*`); `tools:definitions`,
    `tools:skill`, `tools:entry` (key: the group's key; `tools/` `LABEL_*`).
- `Json`, `JsonObject`, `Reference` are in `values.ts`.

## The port (`src/messages/models.ts`)

```ts
interface Provider {
  capabilities(configuration: ModelConfiguration): Promise<ModelCapabilities>;
  project(input: ProjectInput): Promise<ProjectedRequest>; // pure: same input, same payload and digest
  invoke(input: { request: ProjectedRequest; signal: AbortSignal }): AsyncIterable<ModelFrame>;
}
interface ProviderBindings { resolve(target: ModelTarget): Promise<Provider> }
```

- `ModelConfiguration { target: { provider, model }, parameters: { temperature?, topP?, maxOutputTokens?, reasoning? }, cacheScopeId? }`.
  `ReasoningParameters`: `{ mode: 'off' }`, `{ mode: 'adaptive', effort? }`, `{ mode: 'budget', tokens }`. Parameters
  are requests: a provider rejects what its model cannot do.
- `ModelCapabilities`: `contextWindowTokens`, `maxInputTokens?`, `maxOutputTokens?`, `defaultMaxOutputTokens?`,
  `outputLimitIncludesReasoning`, `appendsTools`, `cacheScope`.
- `ProjectInput { messages, tools (available: enabled and loaded), configuration, capabilities }` →
  `ProjectedRequest { payload, digest, revision, inputTokens (an estimate, for admission), prefix? }`. `prefix` is the
  short digests of tools, system, then each message: a request can reuse a cache when the previous one's `prefix` is a
  prefix of its own (the trace's "prefix break"). (`inputTokens` is the providers' calibrated estimate,
  `providers/common/estimate.ts`; the plain four-characters-a-token one is `estimatedTokens` of `messages/`.)
- `ModelFrame`: `text-delta`, `reasoning-delta`, `usage` (totals so far), then exactly one
  `result { message, finish, usage }`, last. `ModelFinish`: `complete`, `tool-calls`, `output-limit`, `filtered`.
- `ProviderUsage`: `inputTokens` (total, cache included), `outputTokens`, `cacheReadInputTokens`,
  `cacheWriteInputTokens`, `reasoningOutputTokens`, `providerRequestId`. The harness converts to `Usage` (uncached,
  cache read, cache write, output, reasoning).
- Failures: throw `ModelError(message, { transient, effect: 'none' | 'unknown', retryAfterMs? })` (defaults: not
  transient, effect `unknown`; pass `effect: 'none'` when the request surely did nothing, so a retry is safe;
  `retryAfterMs` when the provider said how long to wait, which the harness keeps on the attempt's failure);
  `UnrepresentableInput` (code `unrepresentable`) when the provider cannot place something (a system message
  mid-conversation, a binary in a failed tool result, a medium it does not take).

## The providers (`src/providers/`)

A thin layer over the Vercel AI SDK (`ai`, `@ai-sdk/*`; [.v6/ai/README.md](../../../.v6/ai/README.md) decision 7). Three
factories, each its own entry so a host imports only the SDKs it uses
(`@molstar/ai/providers/{anthropic,openai,google}`; all three from `@molstar/ai/providers`):

| Factory                   | Options (besides `resolveApiKey: () => string \| Promise<string>`, `fetch?`, `tokenEstimate?`, `modelLimits?`)            |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `createAnthropicProvider` | `headers?`, `cache?` (default true: cache breakpoints on system, top level and the last immediate tool)                   |
| `createOpenAIProvider`    | `api?: 'responses' \| 'chat'` (default responses), `baseURL?`, `provider?` (default `'openai'`), `label?`, `body?` (chat) |
| `createGoogleProvider`    | none                                                                                                                      |

- **The payload is `streamText`'s own arguments**, as JSON:
  `{ provider, revision, model, toolIds, system?, messages: ModelMessage[], allowSystemInMessages?, tools?: { [wireName]: { description, inputSchema, providerOptions? } }, providerOptions?, maxOutputTokens, temperature?, topP?, appended? }`.
  `project` builds it once (`common/messages.ts`: one `Message` → `ModelMessage` converter, tool calls carrying their
  replay `providerOptions`), digests it with the synchronous `digest`, and gives its `prefix`; `invoke` checks it is
  this provider's (revision, digest, provider: a `ModelError` with effect `none` otherwise), strips the identity fields,
  wraps each tool's schema in `jsonSchema` and passes the rest to `streamText` as it is (`common/stream.ts`). `toolIds`
  maps a returned call's wire name back to its tool id (a call to an unknown name keeps it).
- **Shared core** `common/`: `createProvider(spec: ProviderSpec)` (`provider.ts`: projection, the input estimate,
  segments, invocation), the converter and the context order (`messages.ts`), the stream mapper and the error mapping
  (`stream.ts`), attachments (`attachments.ts`), capabilities (`resolveProviderCapabilities`: `modelLimits[model]`, else
  the catalog, else an error asking for an override). A provider (`anthropic/provider.ts`, `openai/provider.ts`,
  `google/provider.ts`) gives only what differs: the SDK model factory, its parameter checks, reasoning and caching as
  `providerOptions`, its attachment rule, how a call's reasoning is kept and replayed, and (Anthropic, OpenAI) the
  appended layout (`ProviderSpec.native`: which models, and their `project`; OpenAI's in `openai/native.ts`).
- **Any endpoint of a dialect**: `createOpenAIProvider({ api: 'chat', baseURL, provider: 'lab', modelLimits })` serves a
  compatible endpoint under its own provider id. A foreign `baseURL` uses no catalog: limits come only from
  `modelLimits` (keyed by model id), and tools are never appended. Replay identity is the provider id, never a URL.
- **Errors**: any HTTP error has `effect: 'none'`, and 408, 409, 429 and ≥500 are transient; a failure without a status
  (network, mid-stream) is transient with `effect: 'unknown'`; the SDK's `UnsupportedFunctionalityError` is not
  transient, effect `none`; a response's `retry-after-ms` or `retry-after` (seconds or a date) becomes `retryAfterMs`
  (`retryAfterOf`); the key is redacted from messages; the SDK runs with `maxRetries: 0` (retries are the drive's).
- **Revision**: each provider has a `revision` string (`anthropic-ai-sdk-v4`, `openai-responses-ai-sdk-v4`,
  `openai-chat-ai-sdk-v2`, `google-ai-sdk-v4`), in the request and its payload; `invoke` refuses a request of another
  revision, provider or digest, and the harness refuses a retry that projects another. Raise it whenever projection or
  replay rules change.

## Projection rules

- **Collected** (`appendsTools: false`): every available tool goes in the request's tool list; system text may stand
  only before the first non-system message (else `UnrepresentableInput`). The session therefore sends changed guidance
  later as a user message.
- **Appended** (`appendsTools: true`, a catalog quirk: Claude Sonnet 5.5, GPT-6 Luna): a definition loaded mid
  conversation is placed where it loads, so the prefix stays byte-stable as tools are disclosed. Loads follow the head,
  not the callable set. Anthropic: the head's inventory is the tool list, `deferLoading` but for what loads before the
  conversation, and a later load is a system message with `providerOptions.anthropic.toolChanges` (`tool_addition`).
  OpenAI: a load is a placeholder system message (its tools in `providerOptions.molstar.additionalTools`) which the
  provider's `fetch` rewrites into an `additional_tools` developer item, dropping the tool list (it served only to parse
  calls), because the installed SDK has no such input.
- Reasoning maps per provider (Anthropic thinking budget or adaptive effort; OpenAI `reasoningEffort`: Responses sends
  `off` as effort `none`, chat sends no effort; Google Gemini 3+ thinking level or budget); replayed reasoning rides in
  each call's `providerOptions` (`anthropic.thinking`, `openai.reasoning` and `itemId`, `google.thoughtSignature`,
  `chat.text`) and goes back before the call as reasoning parts.
- Attachments: each provider's rule says which binary media it takes where (`AttachmentRule`; Google takes the most,
  OpenAI chat only images in user messages); text types are sent as text (Anthropic: a text document in user content).
  An unsupported one is `UnrepresentableInput` at projection, so the harness can say so before any request; a corrupt
  one (encoding, length, digest) a plain error. A reference's `digest` is `murmur3-128:<32 hex>`, the messages' own hash
  of the decoded bytes (checked when given).
- `cacheScopeId` reaches the provider only where `cacheScope` (OpenAI Responses: `promptCacheKey`).

## The model catalog (`@molstar/ai/models`, `src/providers/models.ts`)

Data, no dependencies (only a `messages/` type). `MODEL_CATALOG: CatalogEntry[]`, each
`{ provider, id, family?, label?, note?, limits?, reasoning?, quirks?, price? }`.

- A **family** entry matches ids that start with it followed by `-`, `.`, `:`, `/`, `_` or `@` (so `gpt-6` covers
  `gpt-6.1-sol`); a model entry matches itself and its dated snapshots (`-YYYYMMDD`, `-YYYY-MM-DD`). Fields merge:
  families (shorter first), then the model. A new version of a family gets its limits, reasoning and quirks; only a
  listed model has a label and price.
- `modelInfo(provider, id, catalog?)`, `listedModels(provider?)`, `defaultReasoning(provider, id)` (what a host passes
  to the session's `thinking`), `dollarsOf(price, usage)` (dollars per call; the long-prompt factor applies to the whole
  call).
- Providers read only limits and the `appendsTools` quirk. Prices, labels and listing are for the host; the session
  knows no models.

## Hosting

A host makes one provider per provider id it has a key for and gives the session a lookup (Prism's assistant:
`packages/prism/assistant/src/viewer/settings-model.ts`, each provider's SDK imported on its first use):

```ts
const providers: Record<string, Provider | undefined> = { anthropic: createAnthropicProvider({ resolveApiKey: () => keys.anthropic }) };
new AgentSession({ models: { providerFor: (id) => providers[id], default: { provider: 'anthropic', model: 'claude-sonnet-5-5' } }, … });
```

## Tests

`src/providers/_test/`: SSE fixtures and a recording `fetch` (`support.ts`); `port.test.ts` (projection pure, digests, a
retained request of another revision, digest or provider refused, `UnrepresentableInput`, tool definitions never text),
one round trip per provider (`anthropic`, `openai`, `openai-chat`: a foreign endpoint, `google`: text, a tool call,
reasoning replay), `prefix.test.ts` (an appended turn's payload and `prefix` extend the previous ones),
`deferred-tools.test.ts` (the appended layout), `attachments`, `capabilities`, `stream` (failure mapping, retry-after,
key redaction), `labels`, `models` (the catalog). No test reaches the network.
