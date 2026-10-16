# Prism assistant: the specification

`@molstar/prism-assistant`: the AI assistant's UI for Prism, opt-in (an app gets it only by depending on it; neither
prism-ui nor the viewer imports it: `ui/src/_test/opt-in-boundary.test.ts`). The agent itself is the viewer's
(`@molstar/prism-viewer/agent`: tools, skills, policies, `createViewerAgent`; [agent.md](../../viewer/_spec/agent.md)),
on `@molstar/ai`; this package is the person's side of it: the model setup and keys, the command menu's Ask row, the
turn card, the Assistant tab, the Inspector. Turns in the record: [journal.md](../../viewer/_spec/journal.md). The code
is the truth for how.

## Entries

| Module                         | What it is                                                                                                    |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------- |
| `spec.ts` (`./spec`)           | `assistant(options)`: the spec's `agent.assistant`, a lazy loader (types only and one dynamic import)         |
| `viewer.ts` (`./viewer`)       | `createAssistant(host)`: the viewer's `AssistantFactory` (viewer-spec.md §2.4), the chunk `assistant()` loads |
| `index.ts`                     | the views and the model types, for stories and other hosts                                                    |
| `keys.ts`                      | `createKeyStore()`: the person's keys, always scrambled when stored (memory, this tab, this browser)          |
| `inspector.ts` (`./inspector`) | the Inspector's page and views, kept out of the index so the tab and the card do not pull it in               |
| `inspector-access.ts`          | the ways in without the Inspector: `InspectorPageId`, `inspectorCommands`                                     |

## Folders

| Folder       | Holds                                                                                                                                                                                                                                                                                                                                                                                                                          |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `viewer/`    | the wiring to the real agent: `settings-model.ts` (providers, keys, the chosen model, own endpoints), `agent-model.ts` (the turn card from the session's hooks, approvals, Undo turn, the history), `inspector-sources.ts` (the Inspector's data from the session and toolkit), `storage.ts` (the keys and allowed domains in the Storage view), `ui.tsx` (the tab with its Inspect rows, the Inspector page and its commands) |
| `views/`     | the Assistant tab, the setup form, the turn card, the history tab (Prism DOM JSX)                                                                                                                                                                                                                                                                                                                                              |
| `inspector/` | the Inspector's views (trace, tools, limits) over `@molstar/ai/trace`, its page and the Inspect rows (`links.ts`)                                                                                                                                                                                                                                                                                                              |
| `model.ts`   | the view models' contracts (`AgentModel`, `AISettingsModel`, `AgentTurn`, `CommandRun`, `RunSource`, `KeyScope`)                                                                                                                                                                                                                                                                                                               |
| `testing/`   | the Inspector's mock data for its unit tests (the stories run the real assistant over a scripted model: `assistant({ provided })`)                                                                                                                                                                                                                                                                                             |

## Contracts

- **Lazy**: `assistant({ provider?, model?, provided? })` imports nothing until the viewer calls it (the tab opened, the
  Ask row used); the assistant's chunk brings the agent, the AI SDK and zod. A page never passes a key.
- **The setup**: Anthropic, Google, OpenAI and one's own OpenAI-compatible endpoints (chat or responses API; https, or
  http on localhost; an optional key; its models typed or fetched from its `/models`; it must allow CORS), one key each.
  The default model (`claude-sonnet-5-5` unless the host names one) is in effect while its provider has a key; the
  person's choice (provider, model) is kept in localStorage apart from the keys, and the first provider with a key
  stands in while the chosen one has none. The model applies from the next request.
- **A host's own models** (`provided`; a function makes them when the assistant loads): listed first, keyed "from the
  server", in effect with no key of the person's (the stories' scripted model).
- **Keys** are the person's, kept for the `KeyScope` they choose: `reload` (memory), `session` (sessionStorage),
  `browser` (localStorage). Stored keys are always scrambled: a per-browser non-extractable AES-GCM key in IndexedDB
  where WebCrypto and IndexedDB are there (against reading or copying a stored value), else an app key in the source
  (AES-GCM, or XOR without WebCrypto: against a naive read). A stored value is opened by trying the per-browser key,
  then the app key. Nothing guards against code running in the page. A key never reaches the model, the record, a
  session file, a log or a notification; the Storage view lists keys by provider and scope only, with Forget.
- **The origin is shared**: every page on the viewer's origin can use stored keys and read kept records; the app runs on
  an origin of its own (app.md §1), and the Storage view names the origin.
- **Allowed domains** (the agent's Allow <host>, [agent.md](../../viewer/_spec/agent.md) §4): for the session, or kept
  in this browser (the prefs' `domains`); listed in the Storage view with Forget.
- **The Ask row** (prism-ui's command menu, the `CommandMenuAssistant` port: `ready`, `propose`, `ask`, `setup`): under
  the list, always there; it takes what is typed with ⌘↵ (Ctrl+↵), or ↵ when no command matches (it is then the
  default). Not set up, it is quiet ("Set up an AI assistant") and opens the Assistant tab.
- **One turn at a time**: `ask` while a turn runs or waits for an approval is refused in words; Stop ends either (the
  viewer's agent settles a parked turn in the record and the history).
- **The turn card** (the shell's agent slot, bottom left of the canvas): the request and its status, the model's words
  as they stream, a retry note, each call in words (title with arguments) with progress, outcome and duration, reads
  quieter; Stop while it runs. At the end, one summary line, the trace link and Undo turn while the checkpoint is held
  (the newest history step is this turn's), "Undone" after. Esc (not from a field, list, menu or dialog) or ✕ folds it
  to a line in the status column; a click opens it again.
- **Approvals** are answered call by call (Allow, Deny; a call with a blocked host: Deny, Allow once, Allow <host> with
  Remember); the turn goes on once every pending call is answered. A denied call is told to the model.
- **History** (`HistoryItem`): the person's runs one by one, an agent's grouped under its request with its checkpoint; a
  run's `RunSource` says who ran it, so an agent's work is never taken for the person's.
- **The Assistant tab** (rail's foot, by Settings): the model in use, the providers' keys (Add, Change, Remove; a key
  from the server has none) with how long each is kept, own endpoints, and the Inspect rows.
- **The Inspector** is a page over the app (`InspectorPageId`, `assistant.page.inspector`) with Trace, Tools and Limits:
  the Assistant tab's Inspect rows (with counts), the turn card's trace (at its turn) and the commands (Assistant trace
  ⌥T, tools, limits; registered once the assistant has loaded) open it at a focus. In the viewer its data is the
  session's (`agent.inspect`): the trace re-read as the harness changes, the tools as the toolkit changes, the limits
  read-only (the session has no setter).
