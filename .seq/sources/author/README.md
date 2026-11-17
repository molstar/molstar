# The author's words

David Sehnal, 2026-11-17, starting the work (lightly condensed; the decisions are quoted).

## The request

`packages/seq` is a new implementation of `packages/sequences`, based on:

- the previous implementation;
- the principles learned about such languages in meso\* (`.meso`, `apps/meso`, `packages/meso` in molstar-ludo);
- the JBrowse developers' reply to SeqViewSpec;
- `@molstar/core/reactive` and `@molstar/ui` where appropriate: "no more rxjs", reusing the reactive model, state and
  the other primitives, and whatever of `packages/ui` is useful.

`.seq` is bootstrapped from the oz principles, on the current patterns "and perhaps some new ones" from the references.

- "this will be a temporary scaffold during development of packages/seq, which will then be deleted. so no references
  to it from code"
- headers: the standard Mol\* header and the author only (no collaboration line)
- "no mentions of base packages/sequences in the implementation, that will be deleted fully when this is done"
- "we should not use .sexpr file type"
- "we don't need any backwards compatibility for this, this is still dev process for the initial library"
- Storybook: "Seq\* (next)" temporarily
- "the initial scaffold should contain setting up the packages/seq alongside its storybook structure, with docs that
  explain how it will work, what are the building blocks etc; and general simplification, new examples, etc."

## Decisions

- **Format**: "yes, let's keep using JSON with s expression printing for easier inspection."
- **Sources**: "for this, i still want this to have a single protocol. but what the spec should have is the ability to
  define a source in a format which an application provides a "runtime" for .. meaning the streaming protocol should be
  standardized, and we convert other formats onto it" (answering the proposal of a portable genome document with
  bigWig and 2bit by URI).
- **The rest, accepted as recommended** ("others as recommended"):
  - packages `packages/seq/{core,viewer,storybook}`, named `@molstar/seqx-*` while `packages/sequences` exists (its
    names are `@molstar/seq-*`), renamed when it is deleted; alternative not taken: one `@molstar/seq` package with
    subpath entries;
  - the scaffold all in `.seq/`, skills read as files and not linked into the harness, one pointer from
    `packages/seq/AGENTS.md`;
  - header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.` with `@author`
    lines only (the request quoted `2019-2022`, read as the template's shape; new files take the current year);
  - modules: oz's `_core`, architecture, representations, ui, testing, and meso\*'s queries and agent-interface; the
    repository's three test lanes, none of the previous implementation's tiers, allocation guards only for the draw
    and hover loops once they exist;
  - from the JBrowse reply, in the first design: self-describing sources, fetch by region (nothing depends on a
    whole-data manifest), strand on every mapping and pair, a published JSON Schema; deferred: layers across
    sections, the smaller additions.

## Later the same day

- **Maps**: asked whether a map could be computational or compressed rather than every position spelled out; answered
  that the resolved form is already run-length segments and that the spec should write each kind in its shortest form
  (CIGAR, exons, offsets, schemes, computed, or by source). "yes, fix this" ([design/building-blocks.md](../../design/building-blocks.md)).
- **Examples and edge cases**: "we could have few example datasets, but should have more views. i would like to keep the
  richness, the editing experiences etc. basically watch out for "edge cases", so the 2nd implmentation doesn't smooth
  over them" ([design/examples.md](../../design/examples.md), [design/edge-cases.md](../../design/edge-cases.md)).
- **Rendering**: "make sure we keep the render ports, the ability render to svg etc. because in the "review reply" Colin
  was able to render it in Jbrowse" ([design/rendering.md](../../design/rendering.md)).

- **UX foundations**: "i would like to have better visual fidelity and nicer UI in general. so we should write "ux
  foundations" at start too, similar to the ones the packages/ui has .. but target to showing the sequence data, when
  to squircles, etc. … the default behavior should be good too .. consistent logic for selections, toggling, etc."
  ([design/ux-foundations.md](../../design/ux-foundations.md)).
- **Nightingale**: "we could also still keep the reference nightingale rendering of a subset as a sanity check for the
  renderer being general" ([design/rendering.md](../../design/rendering.md)).
- **Structure maps and the map model**: "we need ability to map to structural data .. basically SIFTS, but in abstract
  form"; "computed" is how a map is sourced, not its kind; a "computed" map bound to a runtime at viewer composition;
  "a "compressed version" of the mappings too" ([design/maps.md](../../design/maps.md)).

- **JSON Schema**: "i am not sure we need to publish the JSON schema, at least initially. it could end up being bloated,
  and the main use case will will using the "builder" to create the views" ([design/format.md](../../design/format.md)).

- **Grammar**: "the controls grammar should be configurable" ([design/ux-foundations.md](../../design/ux-foundations.md), section 8).

- **Configuration and extensibility**: "basically the rule is: what we can make configurable, we should"; "i am also
  missing notes on "Extensibility" in the storybook" ([design/extensibility.md](../../design/extensibility.md)).

- **More extension points**: "other extensions are also different formats, node kinds, having custom fields"; and the
  reply's edge cases captured too ([design/extensibility.md](../../design/extensibility.md),
  [sources/seqstar/edge-cases.md](../seqstar/edge-cases.md), section 9).

- **Features**: asked "should we also have a "features" document to catalog what we will be implementing in this
  version?"; the proposal accepted: "ok" ([design/features.md](../../design/features.md)).

- **MessagePack**: "maybe msgpack should be the default? it's interchangeable with JSON. and we already spent a lot of
  time getting the builder in packages/sequences understand typedarrays etc?" ([design/format.md](../../design/format.md)).

- **Scope changes**: genome formats, "we should support at least one as proof of concept"; services and structure files
  "should be easy to do at the "application" layer.. should have example for this in storybook"; the Mol\* link "should
  be a showcase to show that the coordiante system translations between mmCIF label/auth works between whatever system
  a sequence uses (many labs for example have "position indexed" sequences, compute predictions for the letters, and
  want to map them onto strucutres!"; markdown links "is about making anything in the view addressable .. so we should
  have that ... we should strobybooks that showcase the addressing stuff" ([design/features.md](../../design/features.md),
  [design/addressing.md](../../design/addressing.md)).

- **Addresses and ids**: "we should be able to address "rects", and the syntax for the address should be path/(query
  :key) for make this generic. meso uses this schema too. and for identifers: use hierarchical addressig? ipv6 llike? and
  we should have "user-defined refs"" ([design/addressing.md](../../design/addressing.md)).

- **Performance**: "let's add "performance" page to the storybook with our targets? both time and memory wise"
  ([design/performance.md](../../design/performance.md)).

- **Composition and IRs**: "expand on how the "spec" will be composed? the views etc .. similar to how meso composes
  many concepts. and do we need intermediate representation the spec compiles/expands to? i think this should
  streamline many things, including mapping the "streaming ports" etc?" ([design/composition.md](../../design/composition.md)).
- **Performance framing**: "the current packages/sequences frames this differently: what is rendered depends on the
  pixels on the screen, not the amount of data ... so our data structures should be chosen appropriately"
  ([design/performance.md](../../design/performance.md)).

- **Layout and ports**: "should "layouting/spacing/etc" be extensible/configurable?"; archives "we should support out
  of the box .. but via a "port" .. in general .. a lot of the functionality should be wired up via "ports""
  ([design/ports.md](../../design/ports.md), [design/ux-foundations.md](../../design/ux-foundations.md)).

- **Testing**: "also what is our testing strategy?" ([design/testing.md](../../design/testing.md)).

- **Golden answers**: "let's take them as guidelines. i think the testing in the old package got really slow and messed
  up the dev process a lot. we should keep fast test in mind too to implement this efficiently"
  ([design/testing.md](../../design/testing.md)).
- **Imports**: "we can import from other mol* parts now, but very carefully and seletively"
  ([design/packages.md](../../design/packages.md)).
- **The tree core**: "yes i think we can bring this to packages/core here. standard headers, since both are mol*"
  ([design/tree-core.md](../../design/tree-core.md)).

- **The tree core's name and units**: `@molstar/core/repr` "ok"; units "ok let's wait"
  ([design/tree-core.md](../../design/tree-core.md)).

- **Precedents and the base library**: "in general, should we get a catalog of precetents to improve on our base
  library?"; then "yes" to building it ([sources/precedents](../precedents/README.md),
  [design/base-library.md](../../design/base-library.md)). Storybook's design pages moved under "Design".

- **Assets**: "reuse asset manger from packages/core/src/util/assets.ts?" ([design/ports.md](../../design/ports.md)).

- **Viewer composition**: "let's be explicit in designs about how the "viewer" will compose? it should be similar to how
  existing "plugins" already composes, and it should configure all the ports, asset mangers, etc"
  ([design/viewer-composition.md](../../design/viewer-composition.md)).

- **After the independent review** (an advisor's read, 2026-11-17): "one" version; the tree core "codec indepedent and
  about abstract patterns"; "this is not just about scope, but about reconciling the architecture, fixing things from
  the "reply", etc"; "rest as recommended" ([design/version-1.md](../../design/version-1.md)). The review also found that
  three of the author's questions had been recorded as decisions (MessagePack as the default, reusing the asset
  manager, configurable layout); the first is now settled the other way (JSON by default, MessagePack an option, as
  recommended), the other two accepted with the rest.

- **Dynamic nodes**: "wouldn't this be needed to have a good streaming protocol though? i would almost lean to do it
  now, because it might be really hard to do later." ([design/dynamic-nodes.md](../../design/dynamic-nodes.md)).

- **Streaming and patches**: "maybe even have top level "Streaming and patches" .. because in a way, these are related
  features? and be able to track what is undoable etc too?" ([design/changes.md](../../design/changes.md)).

- **Second review**: "let's run one more indepedent review with opus, and look for spec and design gaps too"; then
  "Yes, both in order" (fix the contradictions, then write the missing designs in the review's order).

## Reconciliation

- "JSON with S-expression printing" fits §one-recursive-form (the spelling is the project's choice); meso\* chose
  S-expressions, Seq\* chooses JSON: the manifest's Adaptations.
- "a single protocol, other formats converted onto it" is §executors-by-composition and ports (architecture): the
  document names data in its own format, an application binds a runtime for that format, and every runtime speaks one
  protocol to the viewer ([design/sources.md](../../design/sources.md)).
- "temporary scaffold" sits against §permanent-record (records are kept): they are kept for as long as `.seq/` is;
  what must outlive it moves into `_spec/` first.
