# Principles for building large codebases with agentic systems

High-level rules for complex systems (100k to millions of lines) built by one person or a small team working with agents. The rules are split by concern into modules; each rule that has a code-level instantiation links to a **bootstrap pattern** in its module's `bootstrap-patterns/` ([_core/system.md](_core/system.md#code-patterns)), one folder per pattern with a language-neutral README and a TypeScript implementation. Bootstrap patterns are starting implementations of the foundation (wiring, streams, reactive models, commands, the tree, errors, identity, transactions, tests), to take in where they fit; taken in, one becomes the project's own stable pattern. A project's own patterns, stable or sketch, live in its instance; sketches are hardened through designs (§design-first) and prototypes before they are relied on.

All work on a codebase built on these principles is one loop: a session reads a slice of `main`, works through it, and a merge folds the result back ([_core/recurrence.md](_core/recurrence.md)).

**Adopting it in a project:** point the agent at [bootstrap.md](bootstrap.md); it runs the [adopting procedure](_core/patterns/adopting.md).

This folder is a package a project copies and adapts: how it is built (modules, kinds, the `_` rule, ids, skills, adapting) is [_core/system.md](_core/system.md).

## How to read

Every rule, section and pattern has a stable id, cited as `§<id>` (§one-composition-root); ids say nothing about order. Each module has a `README.md` header (kind, status, scope, sources, depends on). Three skills load them, all in `_skills/`: `principles` (this package's `_core`: how work is run), `principles-router` ([SKILL.md](_skills/principles-router/SKILL.md): which module to read, from the project's manifest) and `principles-self-check` (the four-step check, when work changes the principles or designs). The checks every project wires into its gate are in [_checks/](_checks/README.md). A project adds a fourth skill, its instance skill, named after the repository (`_template/instance-skill/`); it loads first and names the project's manifest, which lists these modules and the project's own.

Suggested reading order: architecture, representations and the core workflow first (they shape everything), then the others as needed. Agents loading only part of the rules should take the architecture rules plus the modules for the area they touch.

## Modules (the catalogue)

The package's modules. A project that adopts it keeps its own manifest in its instance, listing these and its own modules, and the router reads that one.

| Module | Kind | Read before | Covers |
|---|---|---|---|
| [architecture](architecture/README.md) | rule set | packages, ports, wiring, errors, commands, performance, access | layout (lib, product packages, apps, services), provider runtime packages, ports, API surfaces (sync/async/streaming), composition root and DI, error codes, commands, performance, access enforcement with application policies, identity and sessions, one consent model for agents' unrequested effects |
| [representations](representations/README.md) | rule set | state, documents, ids, sync, transactions | declarative state as recursive, streamable trees; committed vs transient state; stable addresses, hashing, collaboration; transactions and app-wide optimistic state |
| [ui](ui/README.md) | rule set | UI code, models, design system, Storybook, views | frontend state and logic (view layer, models, responsiveness), design system and Storybook, views as representations |
| [testing](testing/README.md) | rule set | tests, tiers, conformance suites, CI, evals | test tiers, app-level integration and smoke tests, CI, conformance library |
| [_core](_core/README.md) | rule set | any non-trivial work (the `principles` skill loads it) | stable: how the principles work ([system.md](_core/system.md)) and how work is run ([workflow.md](_core/workflow.md)): agent context and specs, design/plan/implement, iteration, keeping design, spec and code honest, meeting the author where they are. |
| [_template](_template/README.md) | (skeletons) | adding a module | one skeleton per kind, for adding a module |

The sources the rules are drawn from and tested against (application ideas, existing codebases) live in the project's `.sources/` (§sources); they point at the rules, never the other way (§self-contained).

## The principles in one page

1. Common code in `lib/` (with a dependency-free base layer in `lib/core/`), product building blocks in `packages/<product>/<module>`, related packages grouped in family folders (`lib/repr/*`) with flat names (`@scope/repr-core`), thin compositions in `apps/` and `services/`; provider SDKs only in `*-<provider>-runtime` packages; common code depends on **ports**, and each app and service wires them in its own **composition root** (§structure, §composition-root).
2. APIs use the right surface (**sync, async or streaming**); state is exposed as read-only **behaviors** (§api-boundaries).
3. State is **declarative data** (one tree, spelled as the project chooses: S-expressions, JSON, ...), recursive and streamable; runtimes/executors interpret it and are bound by composition (§declarative-representations).
4. **Committed** state is the document; **transient** state is an overlay that ends in one commit or a cancel; entities have stable ids and trees have Merkle hashes (§committed-vs-transient, §identity-and-hashing).
5. Large or external data is **referenced**, not embedded (§bulk-data-not-nodes); derived data stays out of the representation (§derived-data-outside).
6. Everything the app can do is reachable by a **command**, with least privilege and a recorded trail; errors have **stable codes** (§stable-error-codes, §command-reachable). Agents reach the app through the same commands, under the same policies.
7. The view layer is thin: models own state, simple actions are bound closures, the UI stays responsive (§frontend-state).
8. Performance is designed in (no accidental quadratics, arenas, zero-alloc hot loops) and enforced by **guard tests** (§performance, §tests).
9. Tests are fast and tiered, with one gate before merge; a small **conformance library** supplies the checks that define correct implementations (§tests, §conformance-suites).
10. **Design first**, plan with a checklist, loop with a `state.md`; process files are reconciled before the close and then kept, never cited, as the work's archive (§design-plan-implement).
11. Spec is the map, code is the territory: intent lives in designs, and disagreement is a **signal** to triage, not a verdict (§agent-context, §recursive-loop).
12. Enforce rules with tools where possible (one linter-formatter, a narrow fallback for the rest, one line width), and never let tools rewrite agent context (`.principles/`, the instance's `.<instance>/`, `.design/`, `.sources/`, `_spec/`, `.execution/`, `.records/`, `AGENTS.md`) (the upstream package's `one-linter-formatter` (not adopted here), §no-auto-format).
13. Changes travel as **transactions**: data, all or nothing, one round trip; optimism lives in the shared store, never in a component (§optimistic-transactions).
14. The base offers mechanisms; **applications set policy and meaning** (access, statuses, conflict rules, budgets). Agents meet the author where they are and reconcile each application against the base in the instance's sources (§author, §access).

## Changing these rules

The rules change the same way everything else does (§recursive-loop): evidence from real work, a short note on what disagreed, and a human decision for anything that moves architecture or principles. Keep them light; if a rule needs a framework to be followed, the rule is wrong.

## Open items

1. Which second language comes first (Rust/WASM for hot paths and repr-core, Python for data and evals?); the workspace is already multi-language (the upstream package's `one-entry-point` (not adopted here), the upstream package's `language-neutral-contracts` (not adopted here)), the patterns are TypeScript only for now.
2. The concrete id scheme and hash function for entities and subtrees (§identity-and-hashing), and whether the default server-ordered collaboration log is enough.
3. Whether to enforce relative test speed automatically (flag tests slower than N times their subject) or leave it to review (§proportionate-tests).
4. Not yet covered by any rule: release/deploy and rollback, configuration and environments, structured logging and tracing beyond the command record (§recorded-invocations).
5. Whether the gaps queue gets a hard cap (no new lane while more than about 30 items wait for a decision); a reference prototype grew it to about 100 while lanes kept landing, and it was left open (§gaps-queue).
6. **Standing tension: proportional review.** A gate's review can add a skeptic, a second agent that tries to refute each finding before it is fixed. When every finding cites the exact line it contradicts, the fixer re-reading that line may be enough; skipping the skeptic then loosens a gate on one agent's judgment. Proportional (a skeptic for findings that rest on judgment, a re-read for line-level contradictions) or always? Carried, not decided (from an adopting project, 2026-11-14). *Later evidence (the same project, 2026-11-17)*: 43 lanes went through a reviewer and a skeptic, and the merges name one refuted finding in 43 reviews without recording the other verdicts, so the evidence cannot settle it. Proposed: merge messages name the refuted findings, so a later fold-back can.
7. **Standing tension: parallel lanes amplify late decisions.** Lanes in worktrees hold up, but they rest on shared foundations staying put; a schema core or a root name that changes mid-wave costs merges into running lanes and a rename pass. §name-the-roots integrates part of it; lanes cannot wait for every decision, and that part stays (from an adopting project, 2026-11-14). *Later evidence (the same project, 2026-11-17)*: it recurred as status documents and version lines, resolved by hand in about half of 29 merges into running lanes; §commit-often now makes those hot files, and the rest stands.
8. **Standing tension: review or check.** Over one adopting project's round (2026-11-17), about 60 of 157 must- and should-fix review findings were in classes a tool could catch (an invalid state landed, a stale cache key, a guard aimed at the wrong path), and every check added costs lane time at every gate; which classes move from review to a check, and at what gate cost, is not settled. The largest class, behavior wrong in ways no rule names, is what review is for. Carried, not decided.
9. **Standing tension: a rule or a design, for one domain's lesson.** A lesson from one application whose shape is any pipeline's may go into the rules (one compiled reading went into an instance's representations module); one whose shape may be the domain's stays in that application's design until a second domain shows it (one writer per value, kept in the adopter's design). Where the line falls is a judgment, carried until a second application shows either.
