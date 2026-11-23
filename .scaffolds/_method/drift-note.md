# Drift: <what disagrees, in a few words>

<!-- The template for a drift note, `execution/<work-id>/drift/<topic>.md` (§drift-triage): code, a spec and a design
     disagree, and neither side is bent to fit. A note that exists is open. It is a gap: asked at the next gate as a yes
     or no with its recommendation, never parked for the close. Resolving it applies the outcome; at the close the
     folder is empty (each note resolved or a line in the TODO). A lane reports drift; the orchestrator writes the note. -->

**Level**: <code-vs-design | spec-vs-code | design-vs-architecture | architecture-vs-principles>. **Proposed
outcome**: <code-wrong | spec-stale | design-wrong | undecided>.

**What**: <what the design or spec says, with its section; what the code does>.

**Where**: <files, and the spec and design sections>.

**Evidence**: <call sites, a test, a measurement, a finding, a commit>.

**Options**: (a) <the design follows the code>; (b) <the code follows the design>; <(c) …>.

**Recommendation**: <the option, and when it would change>.

<!-- Once answered: **Resolved** (the author, <date>, "<words>" or as recommended): <the outcome, and who applies it
     where>. Above code-vs-design, every outcome edits an agreed document, so it is the author's call. -->
