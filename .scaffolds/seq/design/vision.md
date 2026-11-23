# Vision

## Intent

Seq\* is a language for sequence views and a viewer for it. A **spec** (a tree, stored as JSON by default) says what data there is
(sequences, alignments, annotations, the maps between them) and how to show it (sections of tracks of layers). Every
position lives in a **coordinate system**, and **maps** connect systems: a protein onto its gene's codons across
introns, an alignment's rows onto its columns, a design onto its wild type, a chain onto its structure. Data placed on
one system shows on every system it maps to, and hover, selection and focus follow it there.

It serves three readers equally: people writing specs by hand or through a builder, agents writing and patching specs
(the language is small, documented from its own description and validated with every issue reported), and other
viewers (a documented tree in JSON, a resolved view and one data protocol, so a reader needs no Seq\* code to open a
spec). Most specs are written with the typed builder; JSON is what they are stored and exchanged as.

## Simpler this time

The previous implementation proved the concepts and grew around them. This one keeps the concepts
([sources/seqstar](../sources/seqstar/README.md)) and starts from fewer, sharper building blocks:

- **One description of the language** from which validation, defaults, the builder's types, the docs
  and the agent skill are all derived. No hand-written tables beside it.
- **One data protocol**: a source names its data in the data's own format; a runtime the application provides turns
  that format into the protocol ([sources.md](sources.md)). Nothing depends on a whole-data manifest.
- **Pure passes, then models, then views**: spec → resolved model (pure, cached by hash) → layout and draw lists
  (pure) → a backend. The reactive layer sits on top, in DOM-free models (`@molstar/core/reactive`); views are thin
  (`@molstar/ui`).
- **Fewer surfaces**: one default UI (on `@molstar/ui`, with a React adapter if needed), behaviors as a short list of
  replaceable policies, a renderer interface with a handful of members.
- **Performance where it matters**: the draw and hover loops are allocation-free and guarded; everything else is
  plain code.
- **What can be configurable, is** (the author's rule), and every registry is open to extensions passed when the viewer
  is composed; the built-ins use the same points ([extensibility.md](extensibility.md)).
- **Simpler means fewer mechanisms, not fewer cases**: every edge case and editing experience of the previous
  implementation is kept, changed or dropped by decision ([edge-cases.md](edge-cases.md)).
- **No process in the code**: comments cite the package's `_spec` section they follow, nothing else.

## Not goals (for now)

Backwards compatibility with SeqViewSpec 1.0 documents (the author: "we don't need any backwards compatibility");
adapters to other viewers (Nightingale, JBrowse) beyond what the spec, the resolved view and the protocol make possible; a published JSON Schema (deferred, [format.md](format.md)); layers across
sections (deferred, not ruled out by the coordinate model).
