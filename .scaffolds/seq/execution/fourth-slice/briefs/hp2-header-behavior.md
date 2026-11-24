# Lane HP2: the headers behavior and registry kind

**Task**: build the viewer half of [track-header-parts](../../../design/track-header-parts.md) (Decision: Canvas slots or a host's DOM parts; The menu; Touch; The geometry; Undo, history, export; Consequences `seqx-viewer`): the `headerActions` registry kind and its checks, `BuiltInHeaderEntries`, `HeaderFacts`, icons read into paths, the `headers` behavior (hover and tap reveal, slot and chip clicks, the menu and its port `headerMenu`, `headerMenuOverlay()`, the parts and their port `headerParts`, `headerPartsDom()`, the long press, the secondary click under `contextClick` with Shift left to the browser), the cursor's header part (`cursor.prev` at the first position, up and down across headers), `describe.ts`'s header words, `viewer.headers`, `track.show` and `layer.show` over several addresses, the config items. Entries namespaced `<id>:<name>` for a host's ([extension-openings](../../../design/extension-openings.md) rule). The minimal viewer is not grown: the behavior is opt-in.

**Output**: the viewer changes with the tests the design lists (each built-in through its slot, menu item and command alike; a host's entry added and removed with its undo; the menu's data in Node, the overlay in jsdom; a secondary click once per ask, EC-156; geometry published only on change); `_spec` Input, Behaviors, Ports, Commands, Bindings, Configuration. List the public surface for the author's acceptance.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/hp2-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
