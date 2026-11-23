# Skills and tools: the agent interface {#skills-and-tools}

- **Implementations:** [TypeScript](typescript.md)

Rules §skills, §tool-surface: skills explain a domain and are discoverable (index and search), partly generated, and small; the tool surface is two tools per domain (`<domain>_read`, `<domain>_write`) plus `skills_read`, each taking one expression string, dispatched onto the query languages (§query-dsls) and the command registry (§commands).

## Problem

One tool per command gives an agent a long tool list that changes with every release, costs context on every call and
breaks prompt caches. Hand-written skill text drifts from the running system, large reads dump more than the agent can
use, write forms are guessed, and a refusal that only says "no" leaves the agent retrying blindly.

## Shape

- **The tool surface** (§fixed-tool-shape), generated once at startup from the list of domains (`doc`, `scene`, `mol`,
  `data`, `sem`, `agent` §agent-protocol, `wf` orchestration's `workflows` (not adopted here)): fifteen tools, a few hundred tokens. Each takes one input,
  `expr`, a string.

  | Tool | Accepts |
  |---|---|
  | `skills_read` | `(index)`, `(skill <domain> [section])`, `(search "text" [k])` |
  | `<domain>_read` | queries in any of the domain's languages, `(get <locator>)`, `(help ...)`, `(more <cursor>)`, read commands, `(do f1 f2 ...)`; never mutates |
  | `<domain>_write` | commands, `(dry-run <form>)`, `(tx (step [name] <form>) ... )` with `(use name field)`; subject to policy |

- **One expression, two spellings** (§one-recursive-form, meso's `convenience-spellings` (not adopted here)): an S-expression, read leniently with
  the repairs reported (§parser-repairs-slips), or the JSON shorthand of the same tree:
  `{"doc.rename": {"id": "a", "title": "b"}}` is `(doc.rename (id "a") (title "b"))`; arrays are lists. A form's
  attributes become command arguments; an attribute with several values becomes a list.
- **`ToolSession`**, per agent session, created by the transport adapter: `actor`, `policy` (§least-privilege),
  `cursors`, `signal`, `budget` (items and characters per call, §compact-results), and `approved(form)`, set by the
  transport after the user approved that exact form. Approvals answer a harness yield (agent-harness's `harness-yields` (not adopted here)): the drive loop
  gets them wherever the application decides (a transport confirmation, an `approval` message in the agent domain).
- **Reads** (`domainRead`): one snapshot per call. `help` answers from the running registries; `get` reads by locator
  (§semantic-layer); `more` continues a cursor; `do` runs several forms into `(results ...)`; a read-effect command runs
  through the registry; anything else is resolved as a query (§query-dsls) and returns locators.
- **Writes** (`domainWrite`): a command, a `dry-run`, or a `tx` run as one transaction (§batched-calls, §transactions).
- **Results**: `(items ...)`, or over budget `(page (count n) (items ...first budget...) (more <cursor>))`; with repairs,
  `(normalized ...)` and `(repairs ...)` are added.
- **Refusals** (§forgiving-honest-edge): code, message, `retryable`, the refusing layer, and when the dispatcher knows
  one, a corrected form:
  ```
  (error (code "dorothy.bad_input") (message "block.mark takes a text range, not done") (retryable false) (layer validation)
    (try (task.complete (target (at doc (page "welcome") (block "w3"))) (from "v:3c1a"))))
  (error (code "access.denied") (layer policy) (retryable false) (hint "Olive may not edit this page; say so"))
  ```
- **Outline first, with affordances** (§outline-first-reads, §affordances): a read of a large representation returns an
  outline; each record (`kind`, `at`, `excerpt`, `status`, `children`, `hash`, `can`) carries the commands the requester
  may run on it, prefilled, holes as `_`, so a write copies a form instead of guessing one (meso's `documents` (not adopted here)).
- **Skills** (§skills-tested): per domain, a `SKILL.md` (front matter `name` and one-line `description` for harness
  discovery; narrative, workflows, pitfalls) and `examples.sexpr` (runnable examples). `(skill <domain>)` returns the
  static text plus a short generated header (languages, dialects, command count) and a pointer to `(help ...)`, never the
  full catalogs. Skill overrides in session documents are merged over it (agent-harness's `agent-harness` (not adopted here)).
  ```md
  ---
  name: mol
  description: Molecular structure: select, measure and style chains, residues and atoms.
  ---
  # mol
  Read with `mol_read`, change with `mol_write`. For the live vocabulary send `(help functions)`, `(help commands)`, `(help dialects)`.

  ## Workflows
  - Find, inspect, then act: query with `(atoms ...)`, check with `(get <locator>)`, change with a command.
  - Prefer `(dry-run ...)` before destructive commands; they fail with `command.confirmation_required` until the user
    approves that exact form.

  ## Pitfalls
  - Residues have several numberings (`mol.auth`, `mol.label`, sequence positions, external databases). Always write locators
    with their space, and convert with `(sem.convert <locator> <space>)` through `sem_read`; never shift numbers by hand.
  - Carry locators, not bare numbers, between calls.
  ```
- **`(help <topic>)`** (§dynamic-skills-follow-system): topics `commands` (name, effect, description), `functions`
  (name, role, doc, over every query language of the domain), `dialects` (id, lowers to), `errors` (the domain's
  codes), `examples`.
- **Search** (§searchable-skills): a `SkillSearch` port, `search(text, k) -> (loc, title, score)...`, over documents
  derived from skill sections, commands, query functions and error codes (`loc` such as `(command mol.set-color)`).
  Default: fuzzy lexical (trigram similarity, the title weighing twice the text). Later: embeddings behind the same port.
- **A session**:
  ```
  → skills_read  (index)
  ← (skills (doc "Documents: notes and their structure") (mol "Molecular structure: chains, residues, atoms") ...)
  → skills_read  (skill mol)
  ← (skill mol (summary "...") (languages mol.select mol.measure (dialect mol.short)) (read-next (help functions) (help examples)))
  → mol_read     (atoms (residue-name HIS) residue
  ← (items (at mol (chain "A") (res 42)) (at mol (chain "A") (res 97)))   (normalized "...") (repairs (closed-parens "1 at end"))
  → mol_write    (mol.set-color (target (at mol (chain "A") (res 42))) (color "red"))
  ← (ok (op "op_01J..."))
  → skills_read  (search "colour residues")
  ← (matches (command mol.set-color 0.71) (skill mol "Workflows" 0.42) (function mol.select residue 0.38))
  ```
- **Errors**: `tool.wrong_tool` (a mutating form sent to `_read`), `tool.unknown_cursor` (rerun the query).

## Invariants

- The tool list never changes shape after startup: names and input schema are fixed; new commands, functions, concepts
  and dialects appear through `(help ...)`, not as new tools.
- `_read` never mutates; a mutating command sent to it is `tool.wrong_tool` with a hint to use `_write`.
- A read call, `(do ...)` included, sees one snapshot.
- Results stay within the session budget; larger results page with a cursor, never a dump or a silent cut.
- Both spellings of a form give the same result; lenient repairs are reported, never silent.
- `(tx ...)` is all or nothing: one commit per document, or none when any step fails. Commands write into the
  transaction's scope and never commit themselves. When steps need approval, the whole transaction yields once, listing
  them.
- Approval comes only from the transport after the user approved that exact form, never from the agent's own input.
- An unknown command on `_write` is `command.not_found` with a hint to `(help commands)`.
- A refusal carries its code, `retryable` and the refusing layer; the agent stops after a non-retryable refusal.
- An outline record's `can` is filtered by policy for the requester; an empty list is an answer.
- `(help ...)` and the search documents are derived from the running registries, never hand-maintained, so skills cannot
  drift.

## Conformance cases

- Tool list → equals the golden list (names and input schema).
- Every catalog command → reachable through `<domain>_read` or `_write`.
- Any read form → document hash unchanged.
- A large result → within budget, paged with a cursor; an expired cursor → `tool.unknown_cursor`.
- The JSON and S-expression spellings of a form → the same result.
- `(tx ...)` with a failing step → no commit; otherwise one commit per document.
- Every example in a skill → runs through the real tools and matches its expected output.
- Every command, function and dialect a skill names → exists.
- Each skill → fits its context budget. To be decided: the budget (e.g. under 1500 tokens).
- Search → finds every command, function and code by its exact name, and by a one-typo variant.

## Hardened by a prototype

- **Outline first, then drill down**: `page.get` stays a tree for clients and is an outline at the agent's tool edge; a
  page over budget returns its headings with sizes, never a silent cut; derived records name their basis (§outline-first-reads).
- **Affordances**: a record carries the exact write form with one hole; an implied effect's form carries a reason hole
  `(implied "<words>")`, checked against the material; for many writes, offer the bulk form and name addressable parts.
- **Refusals** carry `(try <form>)` for form slips (a wrong field, a stale version with the current one), `retryable`,
  and the refusing layer (enforce, consent, policy, validation, command).
- **Whether vs how**: a capability the agent must know exists (keeping notes beyond this conversation) goes in the
  prompt; the skill alone was not loaded, and the agent said it could not persist anything (§discoverable-skills).

## Open

- Command lookup is unscoped (`mol_read` can run `doc.get`) while `(help commands)` lists only `<domain>.` commands;
  which domain owns commands outside the domain namespaces (`state.query`, §query-dsls; `derive.*`,
  meso's `observable-derivations` (not adopted here)) is not decided.
- The JSON shorthand is not the lossless JSON form of a tree (meso's `convenience-spellings` (not adopted here)); whether the edge accepts both is
  open.
- Whether one value of a list-shaped attribute is read as a single value or a list is left to the command's schema.
- `(help examples)` is listed as a topic, but the help dispatch has no `examples` case yet.
