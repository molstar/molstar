# Drift: language-description.md

## 1. `composeLanguage`

**What**: the design says extensions merge into the language by `composeLanguage(base, entries)`, kinds and params
under their namespace, checked for conflicts, with the builder's static types following a statically composed
language. No `composeLanguage` exists. What is built (extension-openings): `SeqLanguage` is one description;
an extension's layer type is a namespaced word on `layer.type` with its params plain beside the built-in ones (`repr`'s
`extensionParam`), keys `<id>:<name>` ride on any node, an unknown namespaced kind is kept with a warning, and each
extension checks its own (`SeqExtension.check`) where it is composed (core spec, Language, Extensions; repr spec,
Validation).

**Where**: `packages/seq/core/src/lang/language.ts`, `packages/seq/core/src/context/extension.ts` (`SeqExtension`). `packages/seq/TODO.md` lists
`composeLanguage` under `repr`.

**Options**: (a) the design follows the code: the later extension design superseded `composeLanguage`; (b) build it
for extension kinds that need typed builders.

**Recommendation**: (a); (b) only when an extension needs kinds of its own with builder types. If (a) is taken, the
TODO's `composeLanguage` entry is dropped.

## 2. An issue's `path` on a dynamic node

**What**: the design says `path` is the node's place in its answer. The overlay gives the node's form as `address` and,
for the language's checks, a `path` inside the node (`['params', …]`); its own resolve issues carry no path.

**Where**: `packages/seq/core/src/resolve/overlay.ts` (`rules`, `report`); the core spec, Overlay, now says so.

**Options**: (a) the design follows the code (the form names the node; a place in the answer is the source's and no
reader uses it); (b) add the answer's place.

**Recommendation**: (a).

**Resolved** (the author, 2026-11-23, as recommended): item 1: language-description.md (Extensions merge) now describes extensions on the one description, no `composeLanguage`; its row points at core › Extensions. Item 2: the dynamic node's `address` and `path` text amended, its row pointing at core › Overlay. Code item for the orchestrator: drop the `composeLanguage` entry under `repr` in `packages/seq/TODO.md`.
