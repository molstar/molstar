# Stage 1: a spec pass (shared by S1A to S1D)

**Why** (`../plan.md`, stage 1): every agreed design's decisions that the code follows live in a package's `_spec`, so
the scaffold can be archived and still never be cited by code. The specs are for agents: a rule and its reason, terse,
no history (the root `AGENTS.md`, "Keep the docs true").

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/<lane>`.

**For each design your lane brief names, in order**:

1. Read the design's decision sections (`## Decision`, `Decided`, `Defaults`, `Configuration`, the Built notes) and the
   spec(s) your brief names. Read the code where a design and the spec might disagree.
2. Sort each decision into one of three kinds:
   - **in the spec already**: note the section;
   - **followed by the code but missing from the spec**: write it into the spec, at the section where it belongs, as
     a rule and its reason ("Why: …"), in the spec's voice. Never cite `.scaffolds/seq/`, a design file or a principle `§` id;
   - **design-only**: history, intent, a rejected alternative, a plan, a tension still open, or something not built.
     Say which, in a few words.
3. Add an `## In the specs` section at the end of the design (before `## Tensions` if there is one), one row a
   decision:

   ```
   ## In the specs

   | Decision | Where |
   |---|---|
   | <short name> | core › <section> |
   | <short name> | design-only: <reason> |
   ```

   `<spec>` is `core`, `viewer`, `ext-rondo` or `repr`. `<section>` is the start of a heading or a bold lead
   (`- **…**`) in that spec, as written there. Group small decisions into one row when they share a section.
4. Remove the design from `PENDING` in `.scaffolds/seq/_checks/designs-in-specs.ts` and run it: it must pass.
5. Commit per design (or per two small ones).

**Drift**: where the code and the design disagree, the code is what is true now: the spec says what the code does,
the row says "design-only: not built as designed (drift note)", and you write `.scaffolds/seq/execution/sixth-slice/drift/<topic>.md`
(what, where, options, a recommendation). Never change code.

**What you may edit**: the specs your brief names, the designs your brief names (their `## In the specs` section, and
a Built note that is wrong), `designs-in-specs.ts`'s `PENDING`, and `drift/`. Nothing else: no code, no hot file
(`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`), never `state.md`, `plan.md` or `loop.md`. Write only in your
worktree and your own temporary space. Never touch the network.

**Checks**: `pnpm format:check` (run `pnpm prettier --write` on the specs you changed), `pnpm lint`, `pnpm test:unit
packages/seq/ packages/core/src/repr` (the config-table tests read the specs), and the scaffold's checks after
`git add .seq packages`: `node .scaffolds/seq/_checks/designs-in-specs.ts`, `node .scaffolds/principles/_checks/paths-check.ts --allow
'^(?!\.scaffolds/|packages/seq/)'`, `node .scaffolds/principles/_checks/merge-gate.ts --package .scaffolds/principles --instance .scaffolds/seq`, and
`git grep -n '\.scaffolds/' -- packages` (no new hit).

**Report** (under 300 words): commits; per design: rows in the spec, rows design-only, drift notes; the spec sections
you added (names); checks with numbers; anything you could not place.
