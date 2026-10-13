# Modes: the viewer and stories

Prism is one app with more than one mode: what the app is for. **The viewer mode** is Mol\*'s own UI (Home, the
Structure tab, the state tree, Mol\*'s commands). **The story mode** builds and plays stories: a story is a document
(scenes described by schemas, MolViewSpec and SeqViewSpec) that the canvas renders, so it stays serializable and
reproducible; its UX is Stories Untold's rather than Mol\*'s (no state tree, different commands). More modes can follow
the same way. The contracts are in [contexts.md](contexts.md) §1, §3–4.

## 1. One mode per app

An app runs one mode, chosen before the app exists; switching starts a new app. There is no switching inside a running
app, so no state, model or registration crosses from one mode to another, and none of the edge cases of a shared plugin
under two UIs arise.

A mode (`PrismMode`, contexts.md §4) brings, with its contributions:

- its tabs, overlays (cards), sections and toolbar items;
- its commands;
- its modules: the viewer's view models over the plugin; the story mode's document and player;
- its canvas parts: the canvas content, what floats over it (hover label, playback), the canvas stack's lead and aside,
  the context bar; its mouse and touch bindings (the Help tab); whether anything is loaded (`empty`);
- what the agent needs from it: its state block and focus (`account`), and how a turn is undone (`checkpointer`).

The **frame** is the same in every mode: the shell (its layout model, the rail, the docks), the status column, Log
(History with the agent), Settings, Help, the theme, fullscreen, the command menu, and the add-ons (the agent's
Assistant tab, Ask port and turn card).

## 2. Packages and builds: bring what you use

Every mode is its own package, opt-in like the agent: the viewer is `@molstar/prism-viewer` (`packages/prism/viewer`),
the story mode `@molstar/prism-stories` (`packages/prism/stories`), and an app gets a mode only by depending on it. A
build can include one mode only (a story reader without the state tree and the structure tools, a viewer without
stories). prism-ui is the common base they build on and never depends on a mode or names one
(`packages/prism/ui/src/_test/opt-in-boundary.test.ts` checks the modes and the agent). The story mode can split into
`@molstar/prism-stories-*` packages (the document model, the editor, the reader) on the same rule.

## 3. Story mode's model

The story is the source of truth; the plugin renders it. A story command edits the document (add a scene, reorder, edit
a caption, change what a scene shows), and the canvas follows by rendering the scene's MolViewSpec (and the sequence
viewer its SeqViewSpec, [../analysis/harness-and-seqstar.md](../analysis/harness-and-seqstar.md) §2); nothing in the
story mode writes plugin state that the document does not hold. That is what keeps a story reproducible: the same
document renders the same scenes anywhere. Stories are kept in the browser (IndexedDB), so starting another mode loses
none of them.

`StoryModel` (`packages/prism/stories/src/model.ts`): `title`, `scenes` (`StoryScene`: id, title, caption, duration),
`current`, `select`, `addScene` (from what the canvas shows), `removeScene`, and `player` (the scenes as the canvas
plays them; a `SnapshotStepperModel`, so the canvas's playback card shows it; it plays the document, not the plugin's
snapshot manager). It is a module of the story mode (contexts.md §3).

## 4. Choosing a mode

- **The splash**, before any app exists: the wordmark and the modes, each with its icon and a line on what it is for
  (`ModeChooser`, `splash`). Choosing one creates the app in it. A host can name the mode instead
  (`<prism-viewer mode="story">`, `createPrismUI({ mode })`), and the app remembers the last choice (contexts.md §8).
- **Home, while nothing is loaded**, lists the modes (`ModeChooser`, `rows`); choosing another starts a new app in it at
  once (nothing is lost).
- **"Switch mode…"** (the command menu) shows the choice over the running app (`ModeChooser`, `welcome`), saying what
  goes: what is loaded now (stories are kept); choosing disposes the app and creates a new one, Cancel goes back.
- **Clear** (Home's Current session, and the "Clear session" command) empties the session: Home asks once in place, then
  everything loaded goes (the plugin's `clear`; saved snapshots and files stay), and Home offers the modes again.

## 5. State (mocks)

- `PrismAppStory` (`ui/src/testing/shell.tsx`) is the app: the splash, then `ShellStory` (the frame) running one mode;
  "Switch mode…" and Home's choice create it again. `mockViewerMode` (`viewer/src/testing/viewer-mode.tsx`) is Mol\*'s
  own UI as a mode, `ViewerShellStory` the frame with it, `mockBlankMode` the frame alone. A mode's canvas parts are
  render functions, rendered in their slots.
- The story mode (`stories/src/testing`, `mockStoryMode`): the Story tab (the story's title, its scenes, Add scene, the
  current scene's caption), the scenes on the canvas's playback card, story commands (Add scene, Play story, Next and
  Previous scene with Alt+arrows, Export story). Shell › Story mode: started in either mode, and the first open.
- Not yet: the story UX itself (scene editing, what a scene shows, assets, the narrative text, publishing), its document
  model over MolViewSpec and SeqViewSpec, and the modes as `PrismMode`s over `PrismContext` (Phase 6).
