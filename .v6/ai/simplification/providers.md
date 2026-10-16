# Providers: model adapters as a thin layer over the AI SDK

## 1. What the adapter layer does today

`@molstar/ai-model-adapters`: 2,624 source lines, 3,417 test lines (+ `_test/retained-wire.json`, 78). It already runs
on the AI SDK (`streamText`, `@ai-sdk/*`), but around it builds its own projection and replay:

| Part                                   | Files (`model-adapters/src/`)                                               | Lines | What it adds over the SDK                                                                 |
| -------------------------------------- | --------------------------------------------------------------------------- | ----: | ----------------------------------------------------------------------------------------- |
| Shared adapter                         | `common/adapter.ts`                                                         |    90 | `project` (pure: payload, digest, revision, token estimate, segments), `invoke` (checks, then streams) |
| Neutral payload                        | `common/payload.ts`                                                         |   371 | `Message` → a provider-neutral JSON payload (`ProjectedPayload`), measured and digested    |
| Stream and replay                      | `common/stream.ts`                                                          |   417 | decodes the payload again into SDK `ModelMessage`s (`sdkMessages`, `StreamPolicy.inlineMessage`), maps stream parts to `ModelFrame`s, replay state, error mapping |
| Native projections                     | `anthropic/payload.ts` (108), `openai/payload.ts` (100)                      |   208 | mid-conversation tool loads for the `appendsTools` models: Anthropic `tool_addition`, OpenAI `additional_tools`, `deferLoading` |
| OpenAI transport                       | `openai/disclosure-transport.ts`                                            |    51 | a `fetch` wrapper that rewrites placeholder developer messages into `additional_tools` items, because "the installed SDK has no additional_tools input type" (`:13`) |
| Per-provider streams                   | `anthropic/stream.ts` (151), `openai/stream.ts` (137), `google/stream.ts` (61) |   349 | provider options, reasoning replay (signatures, encrypted reasoning), finish mapping       |
| Attachments                            | `common/attachments.ts` (207) + three providers (149)                        |   356 | per-provider media types and placements; `UnrepresentableInput` for the rest               |
| Capabilities                           | `common/capabilities.ts` (75) + three providers (137)                        |   212 | limits from `modelLimits`, else the catalog, else an error; parameter checks               |
| Adapters                               | `anthropic/adapter.ts` (47), `openai/{adapter,chat}.ts` (85, 99), `google/adapter.ts` (43) | 274 | options, SDK model factories, OpenAI chat over `@ai-sdk/openai-compatible`         |
| Estimate                               | `common/estimate.ts`                                                        |    31 | four characters a token                                                                   |
| Catalog                                | `model-catalog.ts`                                                          |   252 | data: limits, reasoning, quirks, prices; `modelInfo`, `listedModels`, `defaultReasoning`, `dollarsOf` |

So a request is converted twice: `Message` → `ProjectedPayload` (retained, digested) → SDK `ModelMessage`s → the SDK's
HTTP body; and for two models a third time inside `fetch`.

**What the harness relies on** (`ModelAdapter` port, `agent/src/models.ts`): `capabilitiesFor(configuration)`
(`harness/src/runner/index.ts:153`), `project(input)` returning `payload`, `digest`, `adapterRevision`, `inputTokens`,
`segments?`, and `invoke` yielding `ModelFrame`s. The runner records `adapterRevision`, `digest` and `segments` as a
model call's evidence and re-projects on retry, refusing a different digest (`harness/src/runner/model.ts:82–92`;
conformance E4 `retry-changed-request`). The trace reads `segments` for prefix breaks (`trace/src/read.ts:310–328`).
The session reads capabilities and wraps adapters (`session/src/session.ts:276,425`).

**Callers**: `packages/prism/assistant/src/viewer/settings-model.ts` (three lazy imports at :207,213,218; the OpenAI one
with `baseURL`, `api`, `provider`, `label` for the person's own endpoints; `lazyAdapter` at :302 implements the port),
`viewer/_test/keyless-endpoint.test.ts` (a `vi.mock` of `@molstar/ai-model-adapters/openai/adapter`),
`viewer/_test/settings-model.test.ts`, `packages/prism/viewer/src/agent/_test/live-models.ts`,
`apps/docs/src/content/docs/v6/prism-viewer/assistant.md:163`.

## 2. The provider, proposed

```ts
interface Provider {
  readonly id: string;                                   // 'anthropic', 'openai', or an endpoint's own id
  capabilities(configuration): Promise<ModelCapabilities>; // the catalog, or the host's limits
  project(input: ProjectInput): Promise<ProjectedRequest>; // straight to streamText's arguments, as JSON
  invoke(input): AsyncIterable<ModelFrame>;               // streamText(payload), parts → frames
}
```

- **Project once, to the SDK's own shape.** The retained payload is `{ model, system, messages: ModelMessage[], tools,
  providerOptions, maxOutputTokens, … }`: JSON, so it still digests, and `invoke` passes it to `streamText` as is. The
  neutral payload and its decode go (`common/payload.ts` mostly, the `sdkMessages`/`inlineMessage` half of
  `common/stream.ts`).
- **One shared converter** `Message` → `ModelMessage` (text, references as file parts, tool calls with their
  `providerOptions` for replay, tool results) and one stream mapper with the error mapping (`retryAfterOf`, transient by
  status, key redaction, `maxRetries: 0`), which is the part worth keeping verbatim.
- **Per provider, only what differs**: the SDK model factory, reasoning → `providerOptions`, Anthropic cache
  breakpoints (`cacheControl` in `providerOptions`), the OpenAI `promptCacheKey`, Google's Gemini 3 rules. Inferred size
  40–100 lines each.
- **Keep the revision and the digest** (one string per provider): the harness's evidence and the retry check depend on
  them, and they cost a line. The digest is the synchronous hash (harness.md §7), so `project` need not be async
  for it.
- **Segments**: computed by the shared code from the projected payload (tools, system, each message), as today in
  `common/adapter.ts:43–51`; ~10 lines.
- **Attachments**: keep the "unsupported media is `UnrepresentableInput`" rule, but let the SDK's own
  `UnsupportedFunctionalityError` decide what is unsupported, mapped in the shared stream code, instead of three
  hand-kept tables (inferred: the SDK raises it for unsupported file parts; to verify per provider).
- **Capabilities**: the catalog's limits, overridden by `modelLimits`; the provider-specific parameter checks become
  "the SDK rejects it" plus the catalog's reasoning modes.

Estimated size (inferred): ~900 source lines including the 252-line catalog; ~1,700 fewer than today.

## 3. The decision: tool additions mid-conversation

The `appendsTools` quirk (catalog `model-catalog.ts:152` Claude Sonnet 5.5, `:177` GPT-6 Luna) places a tool loaded mid
conversation where it loads, so the request's prefix stays byte-stable while skills disclose tools
(`packages/ai/_spec/models.md` § Projection rules). It needs the two native projections (208 lines), the fetch transport
(51), `deferLoading`, and the `deferred-tools.test.ts` (247) and `prefix.test.ts` (287) tests.

- **Keep it**: the providers keep a native path for those two models, and OpenAI keeps the fetch rewrite until the SDK
  has an `additional_tools` input. ~300 lines stay.
- **Drop it**: every available tool goes in the request's tool list (the "collected" rule, already the default for every
  other model). A disclosed tool then changes the tools segment, which providers cache first, so the turn after a
  disclosure re-reads the whole prompt uncached (inferred from the segment order tools → system → messages,
  `common/adapter.ts:43`). The viewer's toolkit discloses skills lazily (`packages/prism/viewer/src/agent/toolkit.ts`), so
  this costs money per disclosure on the default model (`claude-sonnet-5-5` in `packages/ai/_spec/session.md:24`).

Open question 4 in the README.

## 4. Tests

| Test (`model-adapters/src/_test/`)            | Lines | Fate                                                                                   |
| --------------------------------------------- | ----: | -------------------------------------------------------------------------------------- |
| `support.ts` (SSE fixtures, recording fetch)  |   393 | keep, shrink to what the kept tests use                                                |
| `stream.test.ts`                              |   283 | keep: failure mapping, retry-after, key redaction                                      |
| `port.test.ts`                                |   277 | keep a third: projection is pure, digest stable, revision refused                      |
| `anthropic`, `openai`, `google`, `openai-chat` |  1,119 | keep one fixture round trip each (text, a tool call, reasoning replay), ~400            |
| `attachments.test.ts`                         |   351 | keep the `UnrepresentableInput` cases, ~100                                            |
| `prefix.test.ts`                              |   287 | keep one: an appended turn's payload starts with the previous one                      |
| `deferred-tools.test.ts`                      |   247 | goes with §3 "drop"; else stays                                                        |
| `capabilities.test.ts`                        |   143 | keep, shorter                                                                          |
| `compatibility.test.ts` + `retained-wire.json` | 132 + 78 | goes: it pins the neutral payload's bytes by hash, which no longer exist          |
| `model-catalog.test.ts`                       |    99 | keep                                                                                   |
| `labels.test.ts`                              |    86 | keep: labels are never sent                                                            |

Estimated ~1,000 of 3,495 lines stay. These tests are cheap today (all model-adapter files together run in ~0.6 s
locally), so this is about code to maintain, not CI time.

## 5. Renames

`createAnthropicModelAdapter` → `createAnthropicProvider` (and OpenAI, Google), `ModelAdapter` → `Provider`,
`ModelBindings` → `Providers`, `AnthropicModelAdapterOptions` → `AnthropicProviderOptions`, the session's
`models.adapterFor(provider)` → `providers.get(id)` or `providerFor(id)` (`packages/prism/viewer/src/spec.ts:12`
`ModelSource`, the assistant's `settings-model.ts:196`). The harness's `ModelTarget.provider` already says "provider".
