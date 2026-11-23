# <Topic>: <the decision in a line>

<!-- The template for an agreed design, `.scaffolds/<name>/design/<topic>.md` (§design-first). No status, no front
     matter: a file in `design/` is agreed; a draft lives in `execution/<work-id>/drafts/` until the author agrees it.
     Edited in place when it changes (git has the history); never links to a plan. A section marked **Provisional** is
     the agreed direction with its details open, and says what would settle them. -->

## Intent

<!-- What the reader or user can do once this is built, and why. For application-level work (a composition of views,
     navigation, an optional analysis), say so: it is built on the public surface; name the library parts it uses and
     what the library had to gain; the library never ships it as its only way. -->

## Precedents

<!-- The tools whose behaviour this borrows, and where it departs (§precedents); for interaction, what a user can do in
     each that this cannot, each built, planned or refused. A design for a new visible part or gesture names one. -->

- <tool or code>: taken, <what>; not taken, <what>.

## Decision

<!-- What is decided, in sections. The `## In the specs` table maps each decision the code follows to its spec. -->

## Defaults

<!-- For a design that adds a layer, a layout or a part (the defaults-designed check reads Consequences for those
     words): what a reader sees with no setting (sizes, labels, colours, legends; empty, crowded, at the widest and the
     closest view, in dark), each with a sketch or a reference picture, reviewed at the gate. -->

**Mock**: <the mock's stories, one per gesture: crude on purpose, fake data, no performance work; played by the author
at the design gate, linked by each lane's brief, deleted by the merge that lands the built view>. **Precedent**: <the
tool it follows here>.

## Configuration

<!-- One row an item a host or a user can change: copied into the package `_spec`'s Configuration table by the lane
     that builds it. A row names an item the code has, or says "not built". -->

| Item | Default | Why | Alternatives |
|---|---|---|---|

## Rejected alternatives

## Consequences

<!-- What changes elsewhere: the designs and specs it amends, the parts it adds, its edge cases, its performance
     representatives with their paths and tests, its matrix cells (a cell built when a test names its gesture). -->

## Tensions

<!-- What pulls two ways, settled or carried as a standing tension (§surface-tensions). -->

## Self-check

<!-- Run when the design is agreed (§self-check): the design's own claim; what it rests on, negated (for an action,
     usually its inverse or its symmetric case); what the negated world explains; what is built or planned, and what is
     carried as a tension above. -->

## In the specs

<!-- Required once the design has a Decision (the designs-in-specs check): one row a decision the code follows, where
     it lives in a package's `_spec` (`<spec> › <section>`, the spec's key from the instance's map), or
     `design-only: <reason>`. The specs never cite this file. -->

| Decision | Where |
|---|---|
| <decision> | <spec> › <section> |
