# @molstar/seqx-viewer

The viewer: it opens a spec through `@molstar/seqx-core`, keeps the interaction state in models, and draws the core's
draw lists through a backend. It imports `@molstar/seqx-core`, `@molstar/core` and `@molstar/ui`.

## 1. Models own the state

- Every piece of viewer state is in a DOM-free model on `@molstar/core/reactive` (`ReactiveModel`, `State`, `derive`,
  `KeyedState` for per-position and per-row state). Models are tested in Node; timers and frames are injected. Why:
  behavior is testable without a browser, and views stay thin.
- The spec is the committed state; a preview (a drag, a hover's emphasis) is a transient overlay that ends in one patch
  or a cancel. Why: one undo step per gesture, and nothing transient leaks into the spec.

## 2. Views and backends

- Views use `@molstar/ui/dom` and its components (tooltip, context menu, legend, hover label, commands and keymap).
  React, if offered, is an adapter over the same models.
- The Canvas2D backend draws the core's draw lists on screen; export uses the core's SVG backend on the same lists. Why:
  one layout, every output the same picture.
- The draw and hover loops allocate nothing once warm, and are guarded.

## 3. Foundations and grammar

- The data's look extends `@molstar/ui`'s foundations: data tokens (`--seq-*`) derive from `--msui-*`; continuous data
  is square, discrete marks are squircles drawn as paths (square when small or cut); emphasis stacks data, highlights,
  selection, hover, cursor, in that order. Why: one theme, one look in every layout and backend.
- One interaction grammar for every layout and target: a gesture acts on the element under the pointer at the
  granularity hit. The grammar is configuration: a table binding gestures (with modifiers and the kind of element hit)
  to actions, each action a command; the application composing the viewer replaces bindings or picks a preset, and a
  user's changes layer over it. By default click selects (or clears the only selected element), Shift+click extends from
  the anchor, Cmd/Ctrl+click toggles, drag on data selects a snapped range, double-click focuses, Esc clears. One
  configuration applies to every view of a viewer; conflicting bindings are reported when it is built; every action
  stays reachable by a command and the keyboard. Selection is a location projected onto every linked system. Hover,
  selection and the viewport are not history; spec edits are, one undo step each. Why: what a gesture does never depends
  on which view it is in, and an application can still choose what gestures do.

## 4. Composition

- A viewer is made from an explicit `SeqViewerSpec`, as a Mol\* plugin from its `PluginSpec`: registry entries
  (registered in order, atomic, counted, undoable), behaviors with params, bindings, config items, port factories with
  ownership (`disposeWithViewer`), the UI's parts. `DefaultSeqViewerSpec()` lists the defaults; nothing registers on
  import. `createSeqViewer(spec, target?)` creates it, headless without a target; `retire()` disposes what it made and
  undoes its registrations. Why: one composition root shows everything a viewer contains and reaches, and a host shares
  its own ports (a Mol\* plugin's asset manager) without the viewer owning them.

## 5. Configuration and extensions

- The viewer configuration (bindings, behaviors and their params, density, tokens, animation, chrome parts, defaults,
  extensions and runtimes present) is set when an application composes the viewer; user settings layer over it. One
  configuration applies to every view of a viewer. Extensions add commands, behaviors, bindings, structure runtimes and
  chrome parts (track header actions, context menu items, hover card sections, panels) through the same descriptors as
  the core's points.

## 6. Behaviors and commands

- What a gesture does is a behavior: a small policy over the models, replaceable. Everything the viewer can do is a
  command (`@molstar/ui/commands`). Why: an application changes the policy, not the viewer.

## Folders

None yet beyond the entry.
