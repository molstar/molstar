# Messages and models (`@molstar/ai-agent`, `@molstar/ai-model-adapters`)

## Messages (`agent/src/messages.ts`)

- `Message` is `SystemMessage | UserMessage | AssistantMessage`, each `{ role, parts }`. Parts:
  - `ModelVisiblePart`: `{ kind: 'text', text }` or `{ kind: 'reference', reference: Reference }` (a URI with
    `mediaType`, `byteLength`, `digest`; adapters inline only `data:` URIs, others become text `[ref <uri>]`).
  - system only: `ToolDefinitionPart { definition, deferred }`; `deferred: false` loads the tool, `true` only declares
    it. Never sent as text.
  - user only: `ToolResultPart { requestId, tool, status, parts }`.
  - assistant only: `ToolCallPart { requestId, tool, arguments, providerOptions? }` (`providerOptions` keeps the
    provider's replay data: thinking signatures, encrypted reasoning).
- `ToolDefinition { tool: { id }, name, description, inputSchema }`; the provider-safe wire name is derived from the id.
- `hashMessage(message)` (SHA-256 over a structural encoding, async) is an entry's id; `digest(json)` addresses blobs,
  segments and heads; `canonicalJson`, `encodeMessage` are the encodings (`agent/src/hash.ts`).
- `Json`, `JsonObject`, `Reference` are in `values.ts`. Messages carry no metadata: anything the harness needs to know
  about an entry lives elsewhere.

## The port (`agent/src/models.ts`)

```ts
interface ModelAdapter {
  capabilitiesFor(configuration: ModelConfiguration): Promise<ModelCapabilities>;
  project(input: ProjectInput): Promise<ProjectedRequest>; // pure: same input, same payload and digest
  invoke(input: { request: ProjectedRequest; signal: AbortSignal }): AsyncIterable<ModelFrame>;
}
interface ModelBindings { resolve(target: ModelTarget): Promise<ModelAdapter> }
```

- `ModelConfiguration { target: { provider, model }, parameters: { temperature?, topP?, maxOutputTokens?, reasoning? }, cacheScopeId? }`.
  `ReasoningParameters`: `{ mode: 'off' }`, `{ mode: 'adaptive', effort? }`, `{ mode: 'budget', tokens }`. Parameters
  are requests: an adapter rejects what its model cannot do.
- `ModelCapabilities`: `contextWindowTokens`, `maxInputTokens?`, `maxOutputTokens?`, `defaultMaxOutputTokens?`,
  `outputLimitIncludesReasoning`, `appendsTools`, `cacheScope`.
- `ProjectInput { messages, tools (available: enabled and loaded), configuration, capabilities }` →
  `ProjectedRequest { payload, digest, adapterRevision, inputTokens (an estimate, for admission), segments? }`.
  `segments` are short digests of tools, system, then each message: a request can reuse a cache when the previous one's
  segments are a prefix of its own (the trace's "prefix break").
- `ModelFrame`: `text-delta`, `reasoning-delta`, `usage` (totals so far), then exactly one
  `result { message, finish, usage }`, last. `ModelFinish`: `complete`, `tool-calls`, `output-limit`, `filtered`.
- `ProviderUsage`: `inputTokens` (total, cache included), `outputTokens`, `cacheReadInputTokens`,
  `cacheWriteInputTokens`, `reasoningOutputTokens`, `providerRequestId`. The harness converts to `Usage` (uncached,
  cache read, cache write, output, reasoning).
- Failures: throw `ModelError(message, { transient, effect: 'none' | 'unknown' })` (defaults: not transient, effect
  `unknown`; pass `effect: 'none'` when the request surely did nothing, so a retry is safe); `UnrepresentableInput`
  (code `unrepresentable`) when the provider cannot place something (a system message mid-conversation, a binary in a
  failed tool result). There is no structured retry-after; the session parses a wait from the message.

## The adapters (`model-adapters/src/`)

Three factories, on the Vercel AI SDK (`ai`, `@ai-sdk/*`):

| Factory                       | Options (besides `resolveApiKey: () => string \| Promise<string>`, `fetch?`, `tokenEstimate?`, `modelLimits?`)            |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `createAnthropicModelAdapter` | `headers?`, `cache?` (default true: cache breakpoints on system, top level and the last immediate tool)                   |
| `createOpenAIModelAdapter`    | `api?: 'responses' \| 'chat'` (default responses), `baseURL?`, `provider?` (default `'openai'`), `label?`, `body?` (chat) |
| `createGoogleModelAdapter`    | none                                                                                                                      |

- **Any endpoint of a dialect**: `createOpenAIModelAdapter({ api: 'chat', baseURL, provider: 'lab', modelLimits })`
  serves a compatible endpoint under its own provider id. A foreign `baseURL` uses no catalog: limits come only from
  `modelLimits` (keyed by model id), and tools are never appended. Replay identity is the provider id, never a URL.
- **Shared core** `common/`: `createProviderAdapter(spec: ProviderSpec)` (`adapter.ts`), the payload builder and token
  measurement (`payload.ts`), streaming, replay and error mapping (`stream.ts`), attachments, capabilities
  (`resolveProviderCapabilities`: `modelLimits[model]`, else the catalog, else an error asking for an override). Each
  provider folder supplies `capabilities.ts`, `attachments.ts`, `stream.ts` and (Anthropic, OpenAI) a native
  `payload.ts`.
- **Errors**: any HTTP error has `effect: 'none'`, and 408, 409, 429 and ≥500 are transient; a failure without a status
  (network, mid-stream) is transient with `effect: 'unknown'`; the key is redacted from messages; the SDK runs with
  `maxRetries: 0` (retries are the drive's).
- **Revision**: each adapter has a `revision` string in its payload; `invoke` refuses a request of another revision,
  provider or digest. Raise it whenever projection or replay rules change.

## Projection rules

- **Collected** (`appendsTools: false`): every available tool goes in the request's tool list; system text may stand
  only before the first non-system message (else `UnrepresentableInput`). The session therefore sends changed guidance
  later as a user message.
- **Appended** (`appendsTools: true`, a catalog quirk: Claude Sonnet 5.5, GPT-6 Luna): a definition loaded mid
  conversation is placed where it loads (Anthropic `tool_addition`, OpenAI `additional_tools` developer items via
  `openai/disclosure-transport.ts`), so the prefix stays byte-stable as tools are disclosed. Loads follow the head, not
  the callable set.
- Reasoning maps per provider (Anthropic thinking budget or adaptive effort; OpenAI `reasoningEffort`: Responses sends
  `off` as effort `none`, chat sends no effort; Google Gemini 3+ thinking level or budget); replayed reasoning rides in
  each call's `providerOptions`.
- Attachments: per-provider media types and placements (Google takes the most; OpenAI chat only images in user
  messages). An unsupported one is `UnrepresentableInput`; a corrupt one (digest, length) a plain error.
- `cacheScopeId` reaches the provider only where `cacheScope` (OpenAI Responses: `promptCacheKey`).

## The model catalog (`@molstar/ai-model-adapters/model-catalog`)

Data, no dependencies (only an `ai-agent` type). `MODEL_CATALOG: CatalogEntry[]`, each
`{ provider, id, family?, label?, note?, limits?, reasoning?, quirks?, price? }`.

- A **family** entry matches ids that start with it followed by `-`, `.`, `:`, `/`, `_` or `@` (so `gpt-6` covers
  `gpt-6.1-sol`); a model entry matches itself and its dated snapshots (`-YYYYMMDD`, `-YYYY-MM-DD`). Fields merge:
  families (shorter first), then the model. A new version of a family gets its limits, reasoning and quirks; only a
  listed model has a label and price.
- `modelInfo(provider, id, catalog?)`, `listedModels(provider?)`, `defaultReasoning(provider, id)` (what a host passes
  to the session's `thinking`), `dollarsOf(price, usage)` (dollars per call; the long-prompt factor applies to the whole
  call).
- Adapters read only limits and the `appendsTools` quirk. Prices, labels and listing are for the host.
- `@molstar/ai-session` has its own `dollarsOf(table, model, usage)` and `creditsOf` over an app `PriceTable`
  ([session.md](session.md)); the session knows no models.

## Hosting

No host in the repository builds adapters yet; the tests are the examples (`model-adapters/src/_test/support.ts`). A
host makes one adapter per provider it has a key for and gives the session a lookup:

```ts
const adapters: Record<string, ModelAdapter | undefined> = { anthropic: createAnthropicModelAdapter({ resolveApiKey: () => keys.anthropic }) };
new AgentSession({ models: { adapterFor: (p) => adapters[p], default: { provider: 'anthropic', model: 'claude-sonnet-5-5' } }, … });
```

## Tests

`model-adapters/src/_test/`: SSE fixtures and a recording `fetch` (`support.ts`); `compatibility.test.ts` replays
recorded requests and checks their HTTP bodies hash exactly (`retained-wire.json`), so any change to wire bytes shows;
`prefix.test.ts` (cache-stable prefixes), `deferred-tools.test.ts`, `attachments`, `capabilities`, `stream` (failure
mapping, key redaction), `model-catalog`, `openai-chat` (a foreign endpoint). No test reaches the network.
