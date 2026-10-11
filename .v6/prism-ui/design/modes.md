# Modes: Mol\* and stories

Prism is one app with more than one mode: what the app is for at the moment. **The viewer mode** is Mol\*'s own UI
(Home, the Structure tab, the state tree, Mol\*'s commands). **Story mode** builds and plays stories: a story is a
document (scenes described by schemas and MolViewSpec) that the canvas renders, so it stays serializable and
reproducible; its UX is Stories Untold's rather than Mol\*'s (no state tree, different commands). More modes can follow
the same way.

## 1. What a mode is

A mode is a set of contributions activated together and deactivated together:

- its tabs, overlays (cards), sections and toolbar items;
- its commands (the command menu and shortcuts run the active mode's);
- its canvas parts: the canvas content, what floats over it (hover label, playback), the canvas stack's lead and aside,
  the context bar;
- its models: the viewer mode's view models over the plugin; the story mode's story document and player.

What stays across modes is the **frame**: the shell (its layout model, the rail, the docks), the status column, Log
(History with the agent), Settings, the theme, fullscreen, the command menu host, and the add-ons (the agent's Assistant
tab, Ask port and turn card work in every mode).

Switching undoes everything the leaving mode registered (every registration returns its undo; the frame records them)
and activates the next. The rail shows a switch while the app has more than one mode, and each mode has a "Switch to …"
command.

## 2. Packages: bring what you use

Every mode is its own package, opt-in like the agent: the viewer mode is `@molstar/prism-viewer`
(`packages/prism/viewer`), story mode `@molstar/prism-stories` (`packages/prism/stories`), and an app gets a mode only
by depending on it. An app can be composed of story mode alone, with no state tree, no structure tools and none of their
code. prism-ui is the common base they build on (the frame, the components, the command menu) and never depends on a
mode or names one (`packages/prism/ui/src/_test/opt-in-boundary.test.ts` checks the modes and the agent). As it grows it
can split into `@molstar/prism-stories-*` packages (the document model, the editor, the reader), on the same rule.

## 3. Story mode's model

The story is the source of truth; the plugin renders it. A story command edits the document (add a scene, reorder, edit
a caption, change what a scene shows), and the canvas follows by rendering the scene's MolViewSpec; nothing in story
mode writes plugin state that the document does not hold. That is what keeps a story reproducible: the same document
renders the same scenes anywhere.

`StoryModel` (`packages/prism/stories/src/model.ts`): `title`, `scenes` (`StoryScene`: id, title, caption, duration),
`current`, `select`, `addScene` (from what the canvas shows), `removeScene`, and `player` (the scenes as the canvas
plays them; it is a `SnapshotStepperModel`, so the canvas's playback card shows it). In `PrismContext`
([shell.md](shell.md) §6) the story document and its session are models of the stories module.

## 4. Choosing a mode

- **The first open** asks which mode to start in: a card over the empty canvas (`ModeChooser`, `welcome`), each mode
  with its icon and a line on what it is for; the panel stays closed until one is chosen. The real app remembers the
  choice and asks only once.
- **Afterwards** the rail's switch changes modes, and so does Home while nothing is loaded (`ModeChooser`, `rows`): an
  empty session is where a different start makes sense.
- **Clear** (Home's Current session, and the "Clear session" command) empties the session: Home asks once in place, then
  everything loaded goes (the plugin's `clear`; saved snapshots and files stay), and Home offers the modes again.

## 5. State (mocks)

- The mock shell is the frame (`ui/src/testing/shell.tsx`, `ShellStory` with the `modes` the story composes);
  `mockViewerMode` (`viewer/src/testing/viewer-mode.tsx`) is Mol\*'s own UI as a mode, and `ViewerShellStory` the frame
  with it; `mockBlankMode` is the frame alone; a mode registers through a recording shell, so leaving it undoes all of
  it. Its canvas parts are render functions, rendered in their slots.
- Story mode (`stories/src/testing`, `mockStoryMode`): the Story tab (the story's title, its scenes, Add scene, the
  current scene's caption), the scenes on the canvas's playback card, story commands (Add scene, Play story, Next and
  Previous scene with Alt+arrows, Export story). Shell › Story mode shows both modes and the switch.
- Not yet: the story UX itself (scene editing, what a scene shows, assets, the narrative text, publishing), its document
  model over MolViewSpec, and the real modes over `PrismContext` (Phase 6: a mode is a named set of `UIContribution`s,
  [shell.md](shell.md) §5).
