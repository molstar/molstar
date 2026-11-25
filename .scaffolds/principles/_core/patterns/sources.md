# Sources {#sources}

Rules §author, §authors-words, §reconcile-dont-correct, §application-vs-base. **A basic guideline, meant to be
iterated on:** a starting shape for keeping what the principles are drawn from, not a fixed process. Adapt it to
the project, and change it when it stops fitting.

## What a source is

A source is something from outside the principles that they are drawn from or tested against: an application idea
in its author's words, an existing codebase, a person's design docs. It lives in `.sources/<name>/`, one folder each.

| Folder | Holds |
|---|---|
| `.sources/` | what came from outside, in its own terms, and how it maps onto the principles |
| `.design/` | what was agreed |
| `.records/` | how the work got there |

A module that draws on a source names it in words in its header's "Sources" field (a pattern collection's source
folder, in its "source" field); the source's folder in `.sources/` links to the module, not the other way
(§self-contained).

## When to add one

- Someone brings an application idea to build on the principles.
- An existing codebase is studied to learn from it, or to reconcile it against the principles.
- A body of material (docs, screenshots, a person's notes) is distilled into a module.

## What it holds

A starting point; add or drop files as the source needs.

- **`README.md`**: the source in the author's words, or a pointer to where it lives (a repository URL and the commit
  or date it was read at, never a local path); who the author is; what kind of source it is (idea, codebase, docs).
- **`reconciliation.md`**: three lists.
  - **Covered by the base**: what the principles already offer, with rule ids.
  - **Application layer**: each concept that stays in the application, with the base primitives it builds on.
  - **Gaps**: what the base would need.

## How it evolves

- **Gaps travel through the usual loop** (§recursive-loop): a gap becomes a design note, and a design may change a
  rule. The reconciliation links to where each gap went.
- **A gap several sources share** is the strongest case for a new base primitive; a gap only one source has usually
  stays in that application.
- **Update the reconciliation** when the principles change under it (a gap closed, a rule moved), or when a
  prototype teaches something about the source. Its git history is the record of how the fit changed.
- **The author's words stay as written.** They are the reference; corrections go in the reconciliation, never into
  the author's text (§reconcile-dont-correct).
