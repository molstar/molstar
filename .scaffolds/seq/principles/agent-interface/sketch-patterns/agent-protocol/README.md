# Agent protocol: plans, claims and messages as a domain {#agent-protocol}

- **Implementations:** [TypeScript](typescript.md)

Rules §agent-communication: agents and people coordinate through the `agent` domain, a representation with links (meso's `graphs` (not adopted here)), a query language (§query-dsls) and commands (§commands), reached through `agent_read`/`agent_write` (§skills-and-tools) and streamed over the wire (meso's `wire-and-storage` (not adopted here)).

## Problem

Without a shared protocol, agents and people working on the same material step on each other: two agents change the
same chain, a question to a person is lost, a crashed agent forgets what it was doing, and coordination lives in free
text no tool can query. A separate coordination system would duplicate storage, permissions and streaming that the
representations already have.

## Shape

The board is a representation:

```
(board (id "b1")
  (actor (id "a.mol")   (kind agent) (role "structure-editor"))
  (actor (id "a.rev")   (kind agent) (role "reviewer"))
  (actor (id "u.anna")  (kind user)  (role "owner"))

  (plan (id "p1") (owner (link "a.mol")) (goal "Color all histidines in chain A red")
    (step (id "s1") (title "Find histidines in chain A") (status done)
          (reads (at mol (chain "A"))) (result (count 7)))
    (step (id "s2") (title "Color them red") (status in-progress) (after (link "s1"))
          (touches (at mol (chain "A") (res-name HIS))))
    (step (id "s3") (title "Ask for review") (status pending) (after (link "s2"))))

  (claim (id "c1") (holder (link "a.mol")) (on (at mol (chain "A"))) (expires (time "2026-11-13T15:05:00Z")))

  (msg (id "m1") (from (link "a.mol")) (to (link "u.anna")) (kind ask) (deadline (time "2026-11-13T16:00:00Z"))
    (body (approve (dry-run (mol.set-color (target (at mol (chain "A") (res-name HIS))) (color "red"))))))
  (msg (id "m2") (from (link "u.anna")) (reply-to (link "m1")) (kind answer) (body (approved true))))
```

Everything the protocol needs is already in place: ids and ops (§committed-transient-identity), links and integrity
(meso's `graphs` (not adopted here)), locators into other domains (§semantic-layer), time as a core concept (§core-concepts-and-embedding),
persistence and streaming (meso's `wire-and-storage` (not adopted here)).

- **Forms**: `actor` (`kind` agent or user, `role`); `plan` (`owner` link, `goal`, `step` list); `step` (`title`,
  `status`, optional `after` link, `reads` and `touches` as locator lists, optional `result`); `claim` (`holder` link,
  `on` locators, `expires` time: an advisory lease, §claims-are-leases); `msg` (`from` link, `to` a list of actor links
  or `"role:x"` / `"topic:x"`, `kind`, optional `reply-to` and `deadline`, `body` a form of any domain).
- **Step statuses**: `pending`, `in-progress`, `done`, `blocked`, `skipped`, `failed`. **Message kinds**: `ask`,
  `answer`, `inform`, `handoff`, `review`, `approval`.
- **Edges**:

  | Edge | Acyclic | When the target is removed |
  |---|---|---|
  | `after` (step order in a plan) | yes | clear |
  | `owner` (plans need their actor) | | reject |
  | `holder` (an actor's claims go with it) | | cascade |
  | `from` (messages keep their sender, for audit) | | reject |
  | `reply-to` | yes | clear |
  | `to` | | clear |

- **Commands** (`agent_write`): `agent.plan.declare` (goal + steps); `agent.step.set` (status, result, touches; revising
  a plan is more of these); `agent.claim` (`on` locators, `ttl` in seconds; renewing is claiming again);
  `agent.release` (a claim id); `agent.send` (a message).
- **Queries** (`agent_read`): the selection language (§query-dsls) with generators `steps`, `plans`, `msgs` and the
  functions `status`, `owner`, `to`, `unanswered` and `touches-overlaps` (staged: resolves the locator once through the
  semantic graph, so overlap is found across spaces):
  ```
  (steps (status in-progress) (touches-overlaps (at mol (chain "A"))))     ; who is about to change chain A
  (steps (status blocked))                                                 ; what is stuck, across all plans
  (msgs (to "u.anna") (unanswered))                                        ; a person's open asks
  (plans (owner "a.mol"))
  ```
- **An inbox is a subscription**: a subscribe frame (meso's `wire-and-storage` (not adopted here)) with
  `(agent (union (msgs (to "<me>") (unanswered)) (msgs (to "role:<my role>") (unanswered))))`, whose results stream
  changes; in-process, a live query (§query-dsls). Nothing else to build.
- **The agent loop**, the same in a test, an eval or production: declare the plan; query and record the first step's
  result; claim the target (on `agent.claimed`, read the holder's plan, wait or ask); set the step in progress with what
  it touches; act through the domain's `_write`; mark it done; release the claim by the id the claim returned; send a
  review message; the review step stays in progress until the answer arrives in the inbox.
- **The checklist in the UI**: a plan renders as a checklist from the same tree, with no separate UI state
  (§frontend-state).
- **Errors**: `agent.claimed` (retryable), `agent.not_owner`, `agent.ask_expired`, `agent.bad_status`.

## Invariants

- Ownership is checked inside the commands: an `agent` actor changes only entities it owns; users may edit any plan
  they own or that is addressed to them. Anything else is `agent.not_owner`. A policy (§commands) sees the command and
  the actor, not the input, so it only restricts which `agent.*` commands a role may run.
- A claim on targets overlapping a live claim of another actor fails with `agent.claimed`, naming that claim, its holder
  and expiry. A claim is live while `expires` is after now (from the `Clock` port); expired claims never conflict.
- Overlap is decided through the semantic graph, so a claim in one space conflicts with a step touching the same
  things in another.
- The edge rules above hold at commit: `after` and `reply-to` are acyclic, a plan's owner and a message's sender cannot
  be removed, an actor's claims go with it.
- An unknown status or an illegal transition is `agent.bad_status`.
- An ask unanswered past its deadline reports `agent.ask_expired`.
- An inbox delivers a new ask once and drops it when answered.
- Plans and claims live in the store, so after a crash an agent reloads them; nothing is held only in the agent.

## Conformance cases

- An agent changes another actor's plan, claim or message → `agent.not_owner`.
- Two actors claim overlapping targets while the first claim is live → the second gets `agent.claimed`; after it
  expires (fake clock) → the second succeeds.
- A claim in `mol` and a step touching the same atoms via `scene` → found as overlapping.
- A new ask addressed to an inbox → delivered once; answered → dropped.
- An unanswered ask past its deadline → `agent.ask_expired`.
- Crash and reload from the store → the agent has its plan and claims.
- To be decided: the legal status transitions, and their `agent.bad_status` cases.

## Open

- Locator slots name a `sem` language with an `at` form, which §core-concepts-and-embedding does not register yet.
- A claim checks and then inserts, which races with a concurrent claim; whether two overlapping claims may both land, or
  the insert is conditional on the checked base (meso's `conditional-writes` (not adopted here)), is not decided.
- `touches-overlaps` needs the semantic graph to resolve a locator to a canonical space and test overlap; neither
  operation is on §semantic-layer's graph yet.
- Step ids (`s1`) must be unique in the document (§stable-entity-ids) across all plans on the board; whether declaring a
  plan scopes them is open.
- The legal status transitions are not listed (only an example: `done` to `pending` without a reason is illegal).
- The board writes times as `(time "<ISO 8601>")`, while the core time form is `(t <value> <unit>)`
  (§core-concepts-and-embedding); which one a `Time` slot accepts is not settled. The `mol` locators here (`res-name`,
  `res`) belong to the plain `mol` space that is open in §semantic-layer.
