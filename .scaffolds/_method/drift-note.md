# Drift: <what disagrees, in a few words>

<!-- The template for a drift note, `execution/<work-id>/drift/<topic>.md` (§drift-triage): code, a spec and a design
     disagree, and neither side is bent to fit. A note that exists is open. It is a gap: asked at the next gate as a yes
     or no with its recommendation, never parked for the close. Resolving it applies the outcome and deletes the note
     (git keeps it); a note still open at the close becomes a followup. There is no status field. A lane reports
     drift; the orchestrator writes the note. -->

**Level**: <code-vs-design | spec-vs-code | design-vs-architecture | architecture-vs-principles>. **Proposed
outcome**: <code-wrong | spec-stale | design-wrong | undecided>.

**What**: <what the design or spec says, with its section; what the code does>.

**Where**: <files, and the spec and design sections>.

**Evidence**: <call sites, a test, a measurement, a finding, a commit>.

**Options**: (a) <the design follows the code>; (b) <the code follows the design>; <(c) …>.

**Recommendation**: <the option, and when it would change>.

<!-- Once answered, the outcome is applied (code-wrong and spec-stale by an agent; above code-vs-design every outcome
     edits an agreed document, so it is the author's call) and the note deleted; the decision, with the author's words,
     goes into the design or `state.md`'s Decided. -->
