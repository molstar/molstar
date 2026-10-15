# Prism assistant: the specification

`@molstar/prism-assistant`: the AI assistant's UI for Prism, opt-in (an app gets it only by depending on it; neither
prism-ui nor the viewer imports it). The agent itself is the viewer's (`@molstar/prism-viewer/agent`: tools, skills,
policies, `createViewerAgent`), on `@molstar/ai-session`; this package is the person's side of it: the model setup and
keys, the turn card, the Assistant tab, the Inspector. Why it is so:
[.v6/prism-ui/spec/agent.md](../../../../.v6/prism-ui/spec/agent.md) §4–§6,
[journal.md](../../../../.v6/prism-ui/spec/journal.md). The code is the truth for how.

## Entries

| Module                   | What it is                                                                                                |
| ------------------------ | --------------------------------------------------------------------------------------------------------- |
| `viewer.ts` (`./viewer`) | `createAssistant(host)`: the viewer's `AssistantFactory` (viewer-spec.md §2.4); the app passes its import |
| `index.ts`               | the views and the model types, for stories and other hosts                                                |
| `keys.ts`                | `createKeyStore()`: the person's keys, always scrambled when stored (memory, this tab, this browser)      |

## Folders

| Folder       | Holds                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `viewer/`    | the wiring to the real agent: `settings-model.ts` (providers, keys, the chosen model, own endpoints), `agent-model.ts` (the turn card from the session's hooks, approvals, Undo turn, the history), `inspector-sources.ts` (the Inspector's data from the session and toolkit), `storage.ts` (the keys in the Storage view), `ui.tsx` (the tab with its Inspect rows, the Inspector page and its commands) |
| `views/`     | the Assistant tab, the setup form, the turn card, the history tab (Prism DOM JSX)                                                                                                                                                                                                                                                                                                                          |
| `inspector/` | the Inspector's views (trace, tools, limits) over `@molstar/ai-trace`, its page and the Inspect rows (`links.ts`)                                                                                                                                                                                                                                                                                          |
| `model.ts`   | the view models' contracts (`AgentModel`, `AISettingsModel`, `AgentTurn`, `CommandRun`)                                                                                                                                                                                                                                                                                                                    |
| `testing/`   | the Inspector's mock data for its unit tests (the stories run the real assistant over a scripted model: `createAssistantWith({ provided })`)                                                                                                                                                                                                                                                               |

## Contracts

- **One turn at a time**: `ask` while a turn runs or waits for an approval is refused in words; Stop ends either (the
  viewer's agent settles a parked turn in the record and the history).
- **Approvals** are answered call by call (Allow, Deny); the turn goes on once every pending call is answered.
- **Keys** never reach the model, the record, a log or a notification; a stored key is sealed (a per-browser AES key in
  IndexedDB where it can be, else an app key) and listed only by provider and scope.
- **The setup's choice** (provider, model) is kept in localStorage apart from the keys and survives a reload; the first
  provider with a key stands in while the chosen one has none.
- **The Inspector** is a page over the app (`InspectorPageId`, `assistant.page.inspector`): the Assistant tab's Inspect
  rows, the turn card's trace and the commands open it at a focus. In the viewer its data is the session's: the trace
  re-read as the harness changes, the tools as the toolkit changes, the limits read-only (the session has no setter).
