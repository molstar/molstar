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
| spec patch | add a track, hide a layer, set a theme, add a bookmark (a ref), edit a numbering | the spec | yes, one step a batch | yes | yes | the nodes it touched |
| snapshot step | enter a story's step | the spec's snapshots, applied over the base | by leaving the step, not by undo | no (entering is navigation) | the snapshots are | `/snapshots/<id>` |
| preview | a drag resizing a track, a theme picked while hovering a menu | a transient overlay over the spec | cancel (Esc); ends in one patch or nothing | its final patch only | no | as the spec's |
| session | the viewport, the selection, the hover, the cursor, open panels, the view shown | the viewer's models | no (selection: a host may opt into a selection history) | no | no (a host may save it beside the spec) | selection and focus as addresses |
| data arrival | a region loaded, a dynamic node arriving, an entry evicted | the data store | no | no | no (sources are, not their data) | by (source, region) and (source, key) |
| configuration | a binding, a preset, a config item | the viewer's configuration | no (settings history is the host's) | no | no | — |
| registration | an extension registered by a behavior | the viewer's registries | its undo function | no | no | — |

- **Undo is the spec's patch history only** (`@molstar/ui/history`'s `HistoryPort`, with `repr/patch`'s inverses):
  one step a batch, a gesture's preview ending in one step. Nothing else enters it.
- **Patches may name streamed things**: a bookmark on a dynamic node (`/data/genes/(node :key T1)`) is a patch like any
  other; undoing it works whether or not the node is held; the ref waits while it is not ([addressing.md](addressing.md)).
- **Streams never become patches**: arriving data never enters the spec or its history (EC-175); "keep this region" is
  a command that writes a patch (an inline copy, a ref), on purpose.
- **Patches can stream too**: a sequence of patches (an agent's edits, a collaborator's, a host's echo) is applied in
  order, each all or nothing; a patch whose base revision no longer matches is refused or rebased by the host, never
  applied blindly (EC-81, EC-93).
- **Every command says its change kind**: a command's description declares whether it writes a patch, changes session
  state or only navigates, so menus can show what undo will undo and agents can tell an edit from a look.

## Consequences

- The features catalog has one section, "Streaming and patches", with an Undo column ([features.md](features.md)).
- Tests check the table: every command's declared kind matches what it changes; undo touches only spec patches; an
  arrival never changes the spec's hash.
