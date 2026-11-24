# Agent protocol: plans, claims and messages as a domain: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Vocabulary

```ts
const LOC: Slot = { kind: "lang", lang: "sem", form: "at" };
export const STATUS = ["pending", "in-progress", "done", "blocked", "skipped", "failed"] as const;
export const MSG_KINDS = ["ask", "answer", "inform", "handoff", "review", "approval"] as const;

export const agentEdges: EdgeKind[] = [
  { id: "after",    acyclic: true, onTargetRemoved: "clear" },    // step ordering inside a plan
  { id: "owner",    onTargetRemoved: "reject" },                  // plans need their actor
  { id: "holder",   onTargetRemoved: "cascade" },                 // an actor's claims go with it
  { id: "from",     onTargetRemoved: "reject" },                  // messages keep their sender (audit, §least-privilege-full-record)
  { id: "reply-to", acyclic: true, onTargetRemoved: "clear" },
  { id: "to",       onTargetRemoved: "clear" },
];

langs.register({ id: "agent", queries: [{ lang: agentQuery, engine: scanEngine(agentQuery) }], forms: [
  { tag: "actor", slots: { kind: { kind: "atom", type: "string" }, role: { kind: "atom", type: "string" } }, doc: "An agent or a user." },
  { tag: "plan",  slots: { owner: { kind: "link", to: "actor", edge: "owner" }, goal: { kind: "atom", type: "string" },
                           step: { kind: "list", of: { kind: "lang", lang: "agent", form: "step" } } }, doc: "A declared plan." },
  { tag: "step",  slots: { title: { kind: "atom", type: "string" }, status: { kind: "atom", type: "string" },
                           after: { kind: "optional", of: { kind: "link", to: "step", edge: "after" } },
                           reads: { kind: "list", of: LOC },
                           touches: { kind: "list", of: LOC },
                           result: { kind: "optional", of: { kind: "any" } } }, doc: "One checklist item." },
  { tag: "claim", slots: { holder: { kind: "link", to: "actor", edge: "holder" }, on: { kind: "list", of: LOC },
                           expires: { kind: "core", concept: "Time" } }, doc: "An advisory lease on targets." },
  { tag: "msg",   slots: { from: { kind: "link", to: "actor", edge: "from" }, to: { kind: "list", of: { kind: "any-of", options: [
                             { kind: "link", to: "actor", edge: "to" }, { kind: "atom", type: "string" } ] } },  // actor, "role:x" or "topic:x"
                           kind: { kind: "atom", type: "string" }, "reply-to": { kind: "optional", of: { kind: "link", to: "msg", edge: "reply-to" } },
                           deadline: { kind: "optional", of: { kind: "core", concept: "Time" } },
                           body: { kind: "any" } }, doc: "A typed message; the body is a form of any domain." },
]});
// LOC: a locator (§semantic-layer) into any domain. Open: 17 registers no "sem" language with an `at` form yet.
```

## Commands (the `_write` side)

```ts
export const agentCommands = (deps: { doc: Doc; clock: Clock }) => [
  cmd("agent.plan.declare", "mutate", "Declare a plan: goal + steps.", (i: { goal: string; steps: StepInput[] }, ctx) =>
    insertPlan(deps.doc, ctx.actor.id, i)),
  cmd("agent.step.set", "mutate", "Update a step's status/result/touches; revising a plan is just more of these.", (i: StepUpdate, ctx) => {
    ownStep(deps.doc, ctx.actor.id, i.step);                                    // throws agent.not_owner
    return setStep(deps.doc, i);
  }),
  cmd("agent.claim", "mutate", "Claim targets for a while; fails if a live claim overlaps.", async (i: { on: Locator[]; ttl: number }, ctx) => {
    const live = overlappingClaims(deps.doc.committed.value, i.on, deps.clock.now());
    const other = live.find((c) => c.holder !== ctx.actor.id);
    if (other) throw agentErrors.error("agent.claimed", `held by ${other.holder} until ${other.expires}`, { details: [sym("claim"), other.id] });
    // Open: check-then-insert races with a concurrent claim; claims are advisory (§claims-are-leases), but whether two overlapping
    // claims may both land, or the insert is conditional on the checked base (meso's `conditional-writes` (not adopted here)), is not decided.
    return upsertClaim(deps.doc, ctx.actor.id, i.on, deps.clock.now() + i.ttl * 1000);   // renewing = claiming again
  }),
  cmd("agent.release", "mutate", "Release own claims.", (i: { claim: string }, ctx) => removeOwnClaim(deps.doc, ctx.actor.id, i.claim)),
  cmd("agent.send", "mutate", "Send a message.", (i: MsgInput, ctx) => insertMsg(deps.doc, ctx.actor.id, i)),
];

export const agentErrors = defineErrors("agent", {
  "agent.claimed":     { description: "Another actor holds a live claim on an overlapping target.", retryable: true },
  "agent.not_owner":   { description: "Only the owner can change this plan, claim or message.", retryable: false },
  "agent.ask_expired": { description: "The ask passed its deadline without an answer.", retryable: false },
  "agent.bad_status":  { description: "Unknown status or an illegal transition (e.g. done -> pending without a reason).", retryable: false },
});
// helpers assumed: cmd (builds a CommandSpec, §commands), insertPlan/setStep/upsertClaim/insertMsg (ops through one Draft, §committed-transient-identity), Clock (§ports-and-adapters),
// ownStep/removeOwnClaim (throw agent.not_owner), overlappingClaims (live = expires > now)
```

Who may change what: see [README.md](README.md#invariants).

## Queries (the `_read` side)

```ts
// agentQuery: the selection language (§query-dsls) plus a few domain functions; overlap goes through the semantic graph (§semantic-layer)
export const agentQuery = withCore(new QueryLang("agent", langs))
  .register("status", ([s]) => (n) => attrOf(n, "status") === nameOf(s))
  .register("owner",  ([a]) => (n, p, ev) => linkTarget(ev.root, p, "owner") === str(a))
  .register("unanswered", () => (n, p, ev) => tag(n) === "msg" && attrOf(n, "kind") === "ask" && !hasReply(ev.root, idOf(n)!, ev))
  .register("touches-overlaps", ([loc]) => ({
    async prepare(ev) {                                                        // staged: resolve once per evaluation (§query-dsls)
      const target = await sem.resolve(toLocator(loc!), ev);                   // to a canonical space (§semantic-layer)
      return (n) => locatorsOf(n, "touches").some((l) => sem.overlapsSync(l, target, ev));
    },
  }));
// helpers assumed: linkTarget, hasReply, locatorsOf, str, the generators steps/plans/msgs and predicate `to` (like 13's atoms),
// sem (the SemanticGraph, §semantic-layer) with resolve and overlapsSync (cached, uses ev.cache): not on 17's SemanticGraph yet
```

Query examples: see [README.md](README.md#shape).

## An inbox is a subscription

```ts
/** Over the wire (meso's `wire-and-storage` (not adopted here)): a subscribe frame with a query; results frames stream changes. Nothing else to build. */
export const inbox = (send: (f: Frame) => void, me: string) =>
  send({ t: "subscribe", stream: `inbox:${me}`, what: parse(`(agent (union (msgs (to "${me}") (unanswered)) (msgs (to "role:${roleOf(me)}") (unanswered))))`) });
// union: the selection combinator (§query-dsls); `or` combines predicates, not selections. helpers assumed: roleOf
```

In-process, the same thing is `observeQuery(board.effective, engine, q, ctx)` (§query-dsls); see [README.md](README.md#shape).

## Agent loop: declare, claim, act, report

```ts
/** How an agent uses the protocol through its tools; the same sequence in a test, an eval or production. */
export async function editWithPlan(t: Tools, target: string) {
  await t.agent_write(`(agent.plan.declare (goal "Color histidines red")
     (steps (step (id "s1") (title "Find")) (step (id "s2") (title "Color") (after (link "s1"))) (step (id "s3") (title "Review") (after (link "s2")))))`);
  const found = await t.mol_read(`(atoms (chain-id "A") (residue-name HIS))`);
  await t.agent_write(`(agent.step.set (step "s1") (status done) (result (count ${countOf(found)})))`);

  const claim = await t.agent_write(`(agent.claim (on ${target}) (ttl 300))`);   // agent.claimed -> read the plan of the holder, wait or ask
  await t.agent_write(`(agent.step.set (step "s2") (status in-progress) (touches ${target}))`);
  await t.mol_write(`(mol.set-color (target ${target}) (color "red"))`);
  await t.agent_write(`(agent.step.set (step "s2") (status done))`);
  await t.agent_write(`(agent.release (claim "${idIn(claim)}"))`);           // the id the claim returned, never a guess

  await t.agent_write(`(agent.send (to "role:reviewer") (kind review) (body (changed ${target})))`);
  await t.agent_write(`(agent.step.set (step "s3") (status in-progress))`);   // done when the review answer arrives in the inbox
}
// helpers assumed: Tools (the <domain>_read/_write functions from 20), countOf, idIn
// Open: step ids ("s1") must be document-unique (§stable-entity-ids) across all plans on the board; whether declare scopes them is open.
```

## The checklist in the UI

```ts
/** A plan renders as a checklist from the same tree; no separate UI state (§frontend-state). */
export const planModel = (board: TreeState, planId: string) =>
  select(board, (root) => stepsOf(root, planId).map((s) => ({ id: idOf(s)!, title: attrOf(s, "title") as string, status: attrOf(s, "status") as Status })));
// helpers assumed: stepsOf, Status
```

## Conformance

```ts
export const agentProtocolSuite = suite<AgentHost>("AgentProtocol",
  kase("an agent cannot change another actor's plan, claim or message (agent.not_owner)", "quick", async (h) => { /* ... */ }),
  kase("overlapping live claims conflict; expired claims do not", "quick", async (h) => { /* fake clock */ }),
  kase("overlap is found across spaces (claim in mol, step touches the same atoms via scene)", "medium", async (h) => { /* ... */ }),
  kase("inbox subscription delivers a new ask once, and drops it when answered", "quick", async (h) => { /* ... */ }),
  kase("unanswered asks past their deadline report agent.ask_expired", "quick", async (h) => { /* ... */ }),
  kase("after a crash an agent reloads its plan and claims from the store", "medium", async (h) => { /* TreeStore replay, 21 */ }),
);
```
