# Packaging: ten packages into `@molstar/ai`

## 1. What exists today

Ten workspace packages under `packages/ai/*` (`pnpm-workspace.yaml:5`), each with its own `package.json`, `tsconfig.json`
(references between them), inventory entry (`scripts/workspace/inventory.json:376–434`), root `tsconfig.json` reference
(`tsconfig.json:188–215`) and Biome block (`biome.json:94–600`, one per package: no Mol\*, Prism or rxjs; no `node:*`
outside components; zod only in tools; each package's allowed workspace imports). All are `private: true`, version
`6.0.0-dev.0`, with the same `exports` shape: `.` plus a `./*` wildcard (and `./testing`, `./fibers`, `./compaction`,
`./zod`, `./mediabunny` where they exist).

The split encodes the layering `agent → harness → {components, tools, trace} → session` with `model-adapters` beside
it (`packages/ai/_spec/README.md` § Packages and layers). Inside one package that layering stays, as folders and
folder rules (§5).

## 2. The target layout and exports

| Folder (`packages/ai/src/`) | From                                                                  | Notes                                                                 |
| --------------------------- | --------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `messages/`                 | `agent/src/{messages,values,hash}.ts`                                 | `models.ts` (the port) moves to `providers/port.ts`                    |
| `harness/`                  | `harness/src/**`, rewritten as the in-memory loop (harness.md); `components/src/harness/` goes | no storage port, no supervisor |
| `harness/_test/`            | the kept conformance cases (`conformance.md`), a short state test (harness.md §6) | plain tests |
| `providers/`                | `model-adapters/src/**` rewritten (`providers.md`); `model-catalog.ts` → `providers/models.ts` |                                  |
| `tools/`                    | `tools/src/**`                                                        | `zod.ts` per open question 7                                          |
| `session/`                  | `session/src/**` (with `compaction/`)                                 | minus `once.ts`, `prices.ts` (dead-code.md)                           |
| `trace/`                    | `trace/src/**` (with `format.ts`)                                     |                                                                       |
| `threads/`                  | `fibers/src/**`, renamed (`threads.md`)                               | the memory rows only                                                  |
| `testing/`                  | `conformance/src/testing/**`, `conformance/src/session.ts` (`createSession`), `conformance/src/cases/helpers.ts` | no vitest import |
| `narration/`                | `narration/src/**`, only if kept (open question 3)                    | otherwise deleted                                                     |

Proposed `package.json` `exports` (explicit; no `./*`):

```json
{
  ".": "./src/index.ts",
  "./harness": "./src/harness/index.ts",
  "./providers": "./src/providers/index.ts",
  "./providers/anthropic": "./src/providers/anthropic.ts",
  "./providers/openai": "./src/providers/openai.ts",
  "./providers/google": "./src/providers/google.ts",
  "./models": "./src/providers/models.ts",
  "./trace": "./src/trace/index.ts",
  "./testing": "./src/testing/index.ts",
  "./_test/*": null
}
```

(each with the usual `molstar-src` / `types` / `import` conditions, as today). Reasons:

- **The root entry must not pull the AI SDK.** Today the assistant imports each provider lazily
  (`packages/prism/assistant/src/viewer/settings-model.ts:207,213,218`) so a viewer without a key never loads an SDK.
  Providers stay off `.`; the per-provider entries keep one SDK per lazy import.
- **`./models` stays dependency-free** (today `@molstar/ai-model-adapters/model-catalog`, imported eagerly by the
  assistant's model list).
- **`./trace` is separate** because only the Inspector needs it; `./trace` absorbs `./format` (57 lines, no
  dependencies).
- **`./testing`** replaces `@molstar/ai-conformance/testing` for the Prism tests and stories that use `scriptedModel`.
- **`./harness`** is for code that drives `createHarness` directly (the tests, compaction-like helpers); the root
  re-exports only the harness types a host names (`TurnConfiguration`, `Limits`, the reader, …).
- Add `"sideEffects": false` (no AI package sets it today, inferred from the ten `package.json`s): the minimal viewer's
  bundle today carries `@molstar/ai-fibers` although `agent: false` (`examples/minimal-prism-viewer/README.md:107`); one
  package must not make that worse.

Dependencies of `@molstar/ai`: `@molstar/core` (only `data/util/hash-functions`, harness.md §7), `ai`, `@ai-sdk/anthropic`, `@ai-sdk/google`, `@ai-sdk/openai`,
`@ai-sdk/openai-compatible`, `zod`, `tslib` (+ `mediabunny` only if narration stays). Dev: `vitest` (root already has it); `fast-check`
leaves with `head.test.ts` and the S1 case (harness.md §6). `fake-indexeddb` leaves the AI
package (it stays in the catalog: prism-viewer and prism-assistant tests use it, `pnpm-workspace.yaml:110`).

## 3. Name collisions on one root entry

Exported today from two packages under one name (from the export scan; checked in source):

| Name                              | Where                                                                                       | Proposal                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| `dollarsOf`, `ModelPrice`         | `model-adapters/src/model-catalog.ts:35,241` and `session/src/prices.ts:12,26`              | delete the session's (dead-code.md §4); the catalog's in `./models` |
| `REQUEST_MARKER`, `FOCUS_MARKER`  | `session/src/compose.ts:19,21` and `trace/src/read.ts:50,52` (copies)                        | one definition in `session/compose.ts`; trace imports it      |
| `estimatedTokens`                 | `session/src/compaction/plan.ts:460` (an edit's size) and `trace/src/read.ts:152` (text)    | rename the compaction one `editTokens`; one text estimate (duplicates.md §2) |
| `ToolSpec`, `ToolOutcome`         | `conformance/src/testing/tools.ts` and `tools/src/`                                         | the testing ones stay under `./testing` only                  |
| `Clock` vs `FiberClock`           | `harness/src/contracts/runner.ts:37`, `fibers/src/domain.ts`                                | one `Clock` (duplicates.md §5)                                |

`EntryId` (harness) and `FiberEntryId` (fibers) do not collide today, but a rename to `ThreadEntryId` must not shorten
to `EntryId` (threads.md §2).

## 4. What breaks outside `packages/ai`

**Imports**: 57 import sites in 35 files (`git grep "from '@molstar/ai-\|import('@molstar/ai-\|vi.mock('@molstar/ai-"`):

| Area                                   | Files                                                                                                                                                                                  |
| -------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `packages/prism/viewer/src/agent/`     | `session.ts` (3), `toolkit.ts` (2), `tools.ts` (2), `params.ts` (2), `eval.ts`; tests `agent.test.ts`, `eval.test.ts`, `params.test.ts`, `live-models.ts` (2)                         |
| `packages/prism/viewer/src/journal/`   | `record.ts` (3, one a dynamic `import('@molstar/ai-components')` at :774), `kept-records.ts`; tests `record.test.ts` (2), `kept-records.test.ts`                                         |
| `packages/prism/viewer/src/` other     | `spec.ts`, `viewer.ts`, `storage/sources.ts`, `ui/history.tsx` (dynamic import at :92), `storage/_test/sources.test.ts` (2)                                                             |
| `packages/prism/assistant/src/inspector/` | `limits.tsx` (2), `model.ts`, `page.tsx`, `parts.tsx` (2), `tools.tsx` (2), `trace.tsx` (2)                                                                                        |
| `packages/prism/assistant/src/viewer/` | `settings-model.ts` (5: `ModelAdapter`, `listedModels`, three lazy provider imports), `agent-model.ts`, `inspector-sources.ts` (2); tests `agent-model`, `inspect-ui`, `inspector-sources`, `keyless-endpoint` (a `vi.mock` of `@molstar/ai-model-adapters/openai/adapter`), `settings-model` |
| `packages/prism/assistant/src/` other  | `testing/mock-inspector.ts` (2), `_stories/assistant.stories.tsx`                                                                                                                      |
| docs                                   | `apps/docs/src/content/docs/v6/prism-viewer/assistant.md:163`                                                                                                                          |

**Plumbing**:

- `packages/prism/viewer/package.json:51–56,71–72` (six deps, two dev) and `packages/prism/assistant/package.json:41–45,53`
  (five deps, one dev) → one `@molstar/ai` each.
- `packages/prism/viewer/tsconfig.json:44–59` (six references) and `packages/prism/assistant/tsconfig.json:13–34` (five)
  → one `../../ai`.
- `tsconfig.json:188–215` (ten references) → one; `scripts/workspace/inventory.json:376–434` (ten entries) → one;
  `pnpm-workspace.yaml:5` `'packages/ai/*'` → `'packages/ai'`; `pnpm-workspace.yaml:99` comment.
- `biome.json:94–600`: ten blocks keyed by package name → per-folder blocks (§5). Each block's "nothing of Mol\*"
  rule excepts `!@molstar/ai-*` (`biome.json:108–109` and in every block); inside one package that exception goes. No
  Biome rule outside these blocks names the AI packages (prism-ui's opt-in boundary is a test, below).
- `packages/prism/ui/src/_test/opt-in-boundary.test.ts:23–32`: the ten names → `@molstar/ai`.
- `vitest.config.ts:17` (`serial = ['packages/ai/**/*.test.ts']`) keeps working on paths; whether the AI tests still
  need the serial project is conformance.md §4.
- Docs naming the packages: `.v6/ai/{README,harness,plan,tools-and-session}.md`, `.v6/ai/analysis/viewer-mode.md`,
  `.v6/prism-ui/{handoff.md,design/agent.md,spec/README.md,spec/agent.md,spec/journal.md,spec/plan.md}`,
  `packages/prism/assistant/_spec/README.md`, `.v6/plans/changelog.md:96,252`, `examples/minimal-prism-viewer/README.md:107`.

## 5. The layering as folder rules

Each Biome block today restricts a package's imports by package name (`"@molstar/ai-harness"` and so on). After the
collapse the same architecture is folder to folder over relative imports:

| Folder        | May import                                   |
| ------------- | -------------------------------------------- |
| `messages/`   | nothing                                      |
| `harness/`    | `messages/`                                  |
| `providers/`  | `messages/`, the AI SDK, zod (the SDK's peer) |
| `tools/`      | `messages/`, `harness/`, zod                 |
| `trace/`      | `messages/`, `harness/`, `session/compose.ts` (the markers, duplicates.md §1) |
| `threads/`    | nothing (or `messages/` for `Json`)          |
| `session/`    | `messages/`, `harness/`, `tools/`            |
| `testing/`    | `messages/`, `harness/`, `providers/port.ts` |

All folders: no Mol\*, Prism, rxjs, `node:*`, except `messages/` importing `@molstar/core/data/util/hash-functions`
(harness.md §7). Expressing "may import" with Biome's `noRestrictedImports` on relative paths works (verified 2026-10-16 with Biome
2.5.15 on a scratch config): a pattern `**/providers/**` (and `**/providers.js`) flags `../providers/models.js` and
`../../providers/models.js` from a restricted folder, and not a folder's own `./values.js`.

## 6. Size

Removed: 10 → 1 `package.json` (460 → ~50 lines), 10 → 1 `tsconfig.json` (173 → ~15), Biome 511 → est. ~150 lines
(inferred: eight folder blocks of ~15 lines), 9 inventory entries, 9 root references. Moved: all 25,921 source and
11,250 test lines; 196 internal import lines in 126 files become relative (`git grep -c "from '@molstar/ai-"` inside
`packages/ai`).
