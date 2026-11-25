# Commands: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Types (common library)

```ts
export interface Schema<T> {
  parse(input: unknown): T;             // throws on bad input; the registry reports command.invalid_input
  describe(): unknown;                  // JSON Schema, for docs/CLI/agent tools (zod, valibot... all fit)
}

export type Effect = "read" | "mutate" | "destructive";

export interface CommandSpec<I, O> {
  name: string;                         // "doc.rename": namespaced by module
  description: string;
  effect: Effect;
  input: Schema<I>;
  output: Schema<O>;
  /** Every change goes through `ctx.tx`, also under `ctx.dryRun`: the scope is an overlay and nothing is real until the
   *  transaction commits, which runTx decides and commands never see (§transactions). So a command needs no dry-run branch.
   *  Nothing outside the documents happens inside a command: an effect on the world (mail, a provider call) is an
   *  outbox entry written through `ctx.tx` and performed after commit by a worker (§all-or-nothing), so a preview, a replica's
   *  local run or a refused transaction never sends anything. `dryRun` only tells the command this call is a preview.
   *  A mutation that would change nothing says so (a coded result or a refusal), never a bare success (§no-op-writes-say-so). */
  run(input: I, ctx: CommandContext): Promise<O>;
}

export interface CommandContext {
  dryRun: boolean;
  idempotencyKey?: string;              // for retries (§safe-to-automate), never a confirmation
  confirmed?: boolean;                  // set only after the user approved this call
  actor: Actor;                         // §recorded-actors, the same type the access port decides on (§one-policy-port); from the verifier, never input
  signal?: AbortSignal;
  tx?: TxScope;                         // the transaction this call writes into (§transactions); mutating commands never commit
}

/** Who acts (§recorded-actors). Lives in a low shared package so commands, access, auth and tx agree without importing each other. */
export interface Actor {
  readonly id: string;
  readonly kind: "person" | "team" | "agent" | "service";
  readonly on?: Actor;                  // acting for: the verified requester, or for an unrequested effect the person whose
                                        // consent setting authorizes it (§requester-authority, §one-consent-model); recorded and checked
}

/** One record per invocation (§recorded-invocations): audit, replay, debugging, basic metrics. */
export interface InvocationRecord {
  name: string; actor: CommandContext["actor"]; inputHash: string;
  code: string;                         // "ok", an error code, "cancelled", or "error" (uncoded bug)
  ms: number;
}

// codes are declared in the module that throws them (§cataloged-codes, §error-codes); a typo in a code does not compile
export const { error: commandError } = defineErrors("command", {
  "command.not_found":             { description: "No command with this name.", retryable: false },
  "command.invalid_input":         { description: "Input does not match the command's schema.", retryable: false },
  "command.forbidden":             { description: "The policy does not allow this command for this actor.", retryable: false },
  "command.confirmation_required": { description: "Destructive command; needs user confirmation (confirmed: true).", retryable: true },
});
```

## Registry and invoker (transport-agnostic)

```ts
export class CommandRegistry {
  private readonly cmds = new Map<string, CommandSpec<unknown, unknown>>();
  constructor(private readonly record: (r: InvocationRecord) => void = () => {}) {}

  register<I, O>(c: CommandSpec<I, O>): void {
    if (this.cmds.has(c.name)) throw new Error(`duplicate command ${c.name}`);   // programmer error (§programmer-errors)
    this.cmds.set(c.name, c);
  }

  find(name: string): CommandSpec<unknown, unknown> | undefined { return this.cmds.get(name); }

  /** Drives docs, CLI help and the agent's (help commands); `_spec/` links to it instead of copying it (§registry-is-spec). */
  catalog() {
    return [...this.cmds.values()].map((c) => ({
      name: c.name, description: c.description, effect: c.effect,
      input: c.input.describe(), output: c.output.describe(),
    }));
  }

  async invoke(name: string, raw: unknown, ctx: CommandContext, policy?: Policy): Promise<unknown> {
    const t0 = performance.now();
    let code = "ok";
    try {
      const c = this.cmds.get(name);
      if (!c) throw commandError("command.not_found", `no command ${name}`, { hint: "list commands with catalog()" });
      let input: unknown;
      try { input = c.input.parse(raw); }                // schema libraries throw their own errors: give them the code
      catch (e) { throw isCoded(e) ? e : commandError("command.invalid_input", `${name}: ${String(e)}`); }
      await policy?.check(c, ctx, input);              // after parsing: who may do it usually depends on the target
      ctx.signal?.throwIfAborted();
      const out = await c.run(input, ctx);
      try { return c.output.parse(out); }
      catch (e) { throw new Error(`${name} returned output that fails its schema`, { cause: e }); }   // a command bug (§programmer-errors)
    } catch (e) {
      code = isCancelled(e) ? "cancelled" : isCoded(e) ? e.code : "error";   // aborts are not failures (§cancellation-codes)
      throw e;
    } finally {
      // hashJson is total (cycles, BigInt, throwing getters become "unhashable"): a throw here would hide the real error
      this.record({ name, actor: ctx.actor, inputHash: hashJson(raw), code, ms: performance.now() - t0 });
    }
  }
}

/** May be async (an access port may need to look something up, §one-policy-port); gets the parsed input, so it can decide per target. */
export interface Policy {
  check(c: CommandSpec<unknown, unknown>, ctx: CommandContext, input: unknown): void | Promise<void>;
}
/** Every policy must agree; the first refusal wins. */
export const allPolicies = (...ps: Policy[]): Policy => ({
  async check(c, ctx, input) { for (const p of ps) await p.check(c, ctx, input); },
});

/** What models, tools and the CLI hold: actor, defaults and policy bound once. */
export type Invoker = (name: string, input: unknown) => Promise<unknown>;
export const bindInvoker = (r: CommandRegistry, ctx: Partial<CommandContext> = {}, policy?: Policy): Invoker =>
  (name, input) => r.invoke(name, input, { dryRun: false, actor: { kind: "person", id: "local" }, ...ctx }, policy);
```

The composition passes a recorder to keep the records (a log port, or a ring buffer for debugging).

## Defining commands (app or module)

```ts
const { error: docError } = defineErrors("doc", {   // a missing entity is the entity's code, not a command code
  "doc.not_found": { description: "No document with this id.", retryable: false },
});

export const renameDoc = (deps: { store: DocumentStore<DocRecord> }): CommandSpec<{ id: string; title: string }, { id: string; title: string }> => ({
  name: "doc.rename",
  description: "Rename a document.",
  effect: "mutate",
  input: schema<{ id: string; title: string }>(),
  output: schema<{ id: string; title: string }>(),
  async run({ id, title }, ctx) {
    const doc = await deps.store.get(id);
    if (!doc) throw docError("doc.not_found", `no doc ${id}`, { hint: "list ids with doc.list" });
    if (!ctx.dryRun) await deps.store.put({ ...doc, title });   // a plain port without transactions honors dryRun itself;
                                                                // a write through ctx.tx needs no branch (§transactions)
    return { id, title };                            // the same answer for a dry run: what would change
  },
});

// packages/product-x/commands/commands.ts: the product's commands, built from its ports (§ports-and-adapters)
export function registryFor(deps: Deps, record?: (r: InvocationRecord) => void): CommandRegistry {
  const r = new CommandRegistry(record);
  r.register(renameDoc({ store: deps.docs() }));
  /* ... one line per command: doc.create, doc.get, doc.list have the same shape ... */
  return r;
}
```

Commands take ports as dependencies, so they run unchanged on the in-memory composition. Point-in-time state is exposed as `read` commands (`doc.get`, `doc.list`); live observation uses streaming queries (§reachable-queries).

### Representation-valued inputs

A string field that holds a representation is parsed leniently at the agent edge, inside the command (§parser-repairs-slips). The command echoes what it understood:

```ts
/** Lenient parse (§representations): { node, repairs, normalized }; throws a repr.* coded error if it cannot be repaired. The parser is
 *  passed in: the commands library sits below the representation modules and cannot import them (§common-vs-product). */
export const reprArg = (parse: LenientParser, raw: string): Lenient => parse(raw);

// inside a command's run():
const { node, repairs, normalized } = reprArg(parseLenient, input.board);
/* ... use node ... */
return { id, normalized, repairs };
```

## Policy: least privilege per actor

```ts
/** An agent acts within the authority of the person it acts for (§requester-authority): whatever the access port lets both the agent
 *  and `actor.on` do on this target, reads included; no separate allow-list. `on` is the verified requester, or for an
 *  unrequested effect (consent `apply`, a watch) the person whose consent setting authorizes it (§one-consent-model). An agent acting
 *  for no one reads only what it may itself and writes nothing. Destructive commands need the user's confirmation
 *  (§least-privilege). `allow` asks the application's access port (§one-policy-port) for one actor. */
export const agentPolicy = (allow: (a: Actor, c: CommandSpec<unknown, unknown>, input: unknown) => Promise<boolean>): Policy => ({
  async check(c, ctx, input) {
    if (ctx.actor.kind !== "agent") return;
    const on = ctx.actor.on;                           // the full Actor as verified, never rebuilt from an id
    if (!on && c.effect !== "read") throw commandError("command.forbidden", `${c.name}: an agent writes only for someone`);
    if (!(await allow(ctx.actor, c, input)) || (on && !(await allow(on, c, input))))
      throw commandError("command.forbidden", `${c.name} is not allowed for this agent or the person it acts for`);
    if (c.effect === "destructive" && ctx.confirmed !== true)
      throw commandError("command.confirmation_required", `${c.name} is destructive`, { hint: "ask the user, then retry with confirmed: true" });
  },
});
```

Secrets never travel through command inputs, outputs or representations; adapters get them from the composition root.

## Agent / CLI tooling layer (thin adapter)

Agents do **not** get one tool per command: that surface would grow with every command and cost context on every turn. They get two tools per domain (`<domain>_read`, `<domain>_write`) that take command forms and dispatch through this registry, and they discover commands with `(help commands)`. Transport, auth, policy and skill text are app-specific.

```ts
// app: one invoker per agent session; after the user approves a destructive call, the tooling layer
// retries that one call through bindInvoker(registry, { ...agentCtx, confirmed: true }, policy).
const agentCtx = { actor: { kind: "agent", id: agent.id, on: session.requester } } as const;   // requester: the verified Actor (§requester-authority)
const agentInvoke = bindInvoker(registry, agentCtx, agentPolicy(allowByAccessPort));
```

Models submit commands as transactions through the replica (`Send` in [§reactive-models](../../../ui/bootstrap-patterns/reactive-models/README.md), [§transactions](../../../representations/bootstrap-patterns/transactions/README.md)) and never reach around them to a store, so UI and agent behavior cannot diverge (§one-path-to-behavior).

## Tests

```ts
// lib/commands/src/conformance.ts (subpath "./conformance"): the catalog suite every registry runs (§conformance). Plain checks
// that throw, no test-runner import, so the library's runtime entry never pulls in the runner (§conformance, §suites-next-to-ports).
export const commandCatalogSuite = suite<CommandRegistry>("CommandCatalog",
  kase("every command is namespaced, described, with a valid effect and both schemas", "quick", (r) => {
    for (const c of r.catalog()) {
      check(/^[a-z][\w-]*(\.[a-z][\w-]*)+$/.test(c.name), `${c.name}: not namespaced`);
      check(!!c.description, `${c.name}: no description`);
      check(["read", "mutate", "destructive"].includes(c.effect), `${c.name}: bad effect`);
      check(c.input !== undefined && c.output !== undefined, `${c.name}: missing schema`);
    }
  }),
);

// commands.quick.test.ts
register(commandCatalogSuite, () => registryFor(local()), { upTo: "quick" });

// doc-commands.medium.test.ts
it("dry-run does not mutate", async () => {
  const deps = local();
  await deps.docs().put({ id: "a", title: "old" });
  await registryFor(deps).invoke("doc.rename", { id: "a", title: "x" }, { dryRun: true, actor: { kind: "service", id: "test" } });
  expect((await deps.docs().get("a"))?.title).toBe("old");
});
```

Helpers assumed: `schema<T>()` (any schema library), `hashJson` (stable, total hash of the raw input), `check` (throws when false), `defineErrors`/`isCoded`/`isCancelled` (§error-codes), `parseLenient`/`Lenient`/`LenientParser` (§representations), `DocumentStore`/`DocRecord`/`Deps`/`local` (§ports-and-adapters), `suite`/`kase`/`register` (§conformance).
