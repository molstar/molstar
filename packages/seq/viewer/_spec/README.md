# @molstar/seqx-viewer

The viewer: it is composed from an explicit spec (registry entries, behaviors, bindings, configuration, ports), keeps
the interaction state in DOM-free models on `@molstar/core/reactive`, and draws the core's draw lists on a Canvas2D. It
imports `@molstar/seqx-core`, `@molstar/core` and `@molstar/ui`.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## The entry

- `ViewerLanguageVersion` is the core's `LanguageVersion`: the viewer draws the language its core reads.
