# Streaming and patches: how state changes

## Intent

The author, 2026-11-17: "maybe even have top level "Streaming and patches" .. because in a way, these are related
features? and be able to track what is undoable etc too?"

Two things change what a viewer shows: edits to the spec (patches) and data arriving (streams). Both are increments
keyed by address and versioned, both flow into the compiled forms the same way, and both must be clear about what is
undone, saved and shared. One model for every kind of change, with a table that says, for each, where it lives and
whether it is undoable, keeps that clear.

## Decision

### One incremental path

- **A patch** changes spec nodes; their hashes change; the resolved model and view recompile what the changed
  addresses reach ([composition.md](composition.md)).
- **An arrival** changes data store entries; their versions change; draw lists rebuild what read them, and the dynamic
  overlay recomposes what its nodes reach ([data-store.md](data-store.md), [dynamic-nodes.md](dynamic-nodes.md)).
- **Both are increments by address with versions**, so the viewer, an export and another renderer apply differences
  the same way whichever caused them.

### Every change has a kind

| Change | Example | Lives in | Undoable | In history | Saved with the spec | Addressable |
|---|---|---|---|---|---|---|
| spec patch | add a track, hide a layer, set a theme, add a bookmark (a ref), edit a numbering | the spec: a new revision (a commit: the root and its parents) | yes, one step a batch | yes (the commit; its patch is derived) | yes | the nodes it touched |
| snapshot step | enter a story's step | the spec's snapshots: stored patches with no base, applied over the current revision | by leaving the step, not by undo | no (entering is navigation) | the snapshots are | `/snapshots/<id>` |
| preview | a drag resizing a track, a theme picked while hovering a menu | a patch over the current revision, applied again when the revision changes ([patches.md](patches.md), Previews) | cancel (Esc); ends in one patch or nothing | its final patch only | no | as the spec's |
| session | the viewport, the selection, the hover, the cursor, open panels, the view shown | the viewer's models | no (selection: a host may opt into a selection history) | no | no (a host may save it beside the spec) | selection and focus as addresses |
| data arrival | a region loaded, a dynamic node arriving, an entry evicted | the data store | no | no | no (sources are, not their data) | by (source, region) and (source, key) |
| configuration | a binding, a preset, a config item | the viewer's configuration | no (settings history is the host's) | no | no | — |
| registration | an extension registered by a behavior | the viewer's registries | its undo function | no | no | — |
| none | copy a link, read the selection as an address | nothing (a command that only reads; lane O's review, so a menu or an agent never takes it for a move) | — | no | no | — |

- **Seq\* owns no undo stack; the host does** (the author, 2026-11-18: the viewer lives inside a bigger app, a
  structure viewer, Prism, a Mol\* plugin, which keeps one undo stack for everything). Seq\* exposes a **history port**
  (`seqx-viewer/ports`) in its own terms, with two faces:
  - *the host drives Seq\** (a participant): `checkpoint()` is a commit, Seq\*'s own type (`SpecCommit`: the spec's
    immutable root, its revision, its parents), cheap by structural sharing; `restore(c)` sets the root back (a move of
    the pointer) and drops a preview; `changed` tells of edits; an optional `root(revision)` lets a host resolve a
    revision older than Seq\*'s recent roots (N-14);
  - *Seq\* asks the host* (a run): an edit runs as a step with its label and change kind, so the host's stack
    interleaves Seq\* edits with its own. One step a batch, a gesture's preview ending in one step.

  Only spec patches enter it. `repr/patch`'s inverses serve `diff`, patch streams and agents, not the undo stack.
  Adapters live in each composition root: `@molstar/ui/history` for the standalone viewer and Prism, the Mol\* plugin's
  state history later, none headless.
- **Patches may name streamed things**: a bookmark on a dynamic node (`/data/genes/(node :key T1)`) is a patch like any
  other; undoing it works whether or not the node is held; the ref waits while it is not ([addressing.md](addressing.md)).
- **Streams never become patches**: arriving data never enters the spec or its history (EC-175); "keep this region" is
  a command that writes a patch (an inline copy, a ref), on purpose.
- **Patches can stream too, in version 1 from one writer at a time**: an agent's or a host's patches are applied
  through the same history, each a step (an agent's turn may group its calls into one, as `HistoryStep.turn` allows);
  a patch names the revision it was made against and is refused when that is not the current one (EC-81), never
  rebased automatically. A host may rebase it on purpose: `merge3(base, current, apply(base, p))`, conflicts shown or
  rescued as data (`merge3`, its own agreement before the editing lane). A change made outside a run is reported to the host through the
  port, and whether it invalidates the earlier steps is the host's rule (`@molstar/ui/history`'s adapter invalidates
  them). Concurrent writers (EC-93) have their primitive in `merge3`; who merges and when is the host's policy, not yet
  designed.
- **Revisions are the truth, patches derived** (git's model, the author 2026-11-18; [patches.md](patches.md)): what a
  history records is a commit, a revision's root with its parents, made by `SpecModel` per applied change and handed to
  the host as the port's checkpoint (N-14); an undo
  restores an earlier root (a move of the pointer, no new commit); the patch from a parent is `diff(parent, root)`
  whenever one is wanted. Why: a recorded root is exact by construction (ref rewrites and normalization included), and
  checkpoints by root cost only the paths that changed. *Snapshot* means a story step only, never a revision.
- **Story steps are stored patches** with no base, in `patchLanguage(SeqLanguage)`, each a full state of the view
  stored as a patch; their refs are typed by kind and key and resolved over the step's dry run ([patches.md](patches.md));
  a step that no longer applies is reported on its node and shown as not applying; edits while a step is shown patch
  the base. How a story's steps are reached, transitioned and animated, with the application owning the clock, is
  [stories.md](stories.md).
- **What is shown with a step and a preview** (N-13): `shown = step(applyPatch(base, preview))`. The preview is an edit
  to the base, and the step is applied over it.
  - An edit or a preview whose target only the step made is refused (`repr.patch.target`, its hint naming the step).
  - One that sets a key the step also sets lands in the base and changes nothing visible: `SpecModel` reports it as
    masked (an info naming the step), never silently.
  - An agent reads the base revision, with the step shown named beside it.
- **Rebase needs the base root**: `SpecModel` keeps the recent roots by revision (`revisions.keepRecent`, default 64),
  and a host may resolve older ones through the port's `root(revision)`.
- **Every command says its change kind**: a command's description declares whether it writes a patch, changes session
  state or only navigates, so menus can show what undo will undo and agents can tell an edit from a look.

## Consequences

- The features catalog has one section, "Streaming and patches", with an Undo column ([features.md](features.md)).
- Tests check the table: every command's declared kind matches what it changes; undo touches only spec patches; an
  arrival never changes the spec's hash.
