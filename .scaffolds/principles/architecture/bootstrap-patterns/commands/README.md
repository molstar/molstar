# Commands {#commands}

- **Implementations:** [TypeScript](typescript.md)

Rules §command-reachable, §stable-error-codes: everything the app does is a named, typed, self-describing command in one registry; UI, tests, CLI and agents all go through it, with a least-privilege policy, coded errors and a record of every invocation.

## Problem

When the UI, the API, the CLI and agent tools each carry their own logic, they drift apart, and an agent operates a
simulation of the app rather than the app. Behavior that is not described as data cannot generate docs, a CLI or a
tool list, cannot be checked by a policy before it runs, and leaves no record of who did what.

## Shape

- **Command spec**: `name` (namespaced by module, `doc.rename`), `description`, `effect` (`read`, `mutate` or
  `destructive`), an `input` and an `output` schema, and `run(input, ctx) -> async output`.
- **Schema**: `parse(raw) -> value` (raises on bad input) and `describe() -> machine-readable schema` (for example JSON Schema).
- **`CommandContext`**: `dryRun`, optional `idempotencyKey` (for retries, never a confirmation), optional `confirmed`
  (set only after the user approved this call), `actor`, an optional cancellation signal, and the optional transaction
  scope `tx` the command writes into (§transactions).
- **`Actor`**: `id`, `kind` (`person`, `team`, `agent`, `service`) and optional `on`, the actor an agent acts for
  (§requester-authority, §one-consent-model); one type shared by commands, access, auth and transactions.
- **Invocation record**: `name`, `actor`, `inputHash`, `code` (`ok`, an error code, `cancelled`, or `error` for an
  uncoded bug), `ms`.
- **`CommandRegistry`**: `register(spec)`, `find(name)`, `catalog()` (name, description, effect and both schema
  descriptions per command; it drives docs, CLI help and `(help commands)`), and `invoke(name, raw, ctx, policy?)`.
  It takes a recorder that receives every invocation record.
- **Policy**: `check(spec, ctx, parsedInput)`, which returns or raises a coded error and may be async; `allPolicies(...)`
  requires every one to pass.
- **`Invoker`**: `(name, input) -> async output`, what models, tools and the CLI hold. `bindInvoker(registry, ctx, policy)`
  binds the context and policy once; the default context is `dryRun: false` and the local person actor.
- **`registryFor(deps)`**: the product's registry, one command per line, built from its ports (§ports-and-adapters).
- **`reprArg(parse, raw)`**: a representation-valued string input parsed leniently inside the command; the command
  echoes what it understood (`normalized`, `repairs`) (§parser-repairs-slips).
- **`agentPolicy(allow)`**: least privilege for agents, asking the access port (§one-policy-port) once per actor.
- **Agent surface**: two tools per domain, `<domain>_read` and `<domain>_write`, that take command forms and dispatch
  through the registry; discovery is `(help commands)`.

Codes this pattern declares (§error-codes):

| Code | Meaning | Retryable |
|---|---|---|
| `command.not_found` | no command with this name | no |
| `command.invalid_input` | input does not match the command's schema | no |
| `command.forbidden` | the policy does not allow this command for this actor | no |
| `command.confirmation_required` | destructive; needs the user's confirmation | yes |

`invoke` runs in this order: look up the name, parse the input, check the policy with the parsed input, check
cancellation, run, check the output against its schema, and record the invocation whatever happened.

## Invariants

- Registering a second command with a name already registered is a programmer error.
- An unknown name raises `command.not_found` with a hint naming how to list commands.
- Input that fails its schema raises `command.invalid_input`; a coded error raised by the schema passes through unchanged.
- The policy sees the parsed input and runs before the command; under `allPolicies` the first refusal wins.
- Output that fails its schema is a programmer error (plain exception), not a coded error.
- Every `invoke` produces exactly one invocation record, also for unknown names and failures; computing the input hash
  never raises.
- A cancelled call is recorded as `cancelled`, not as a failure (§cancellation-codes).
- A dry run returns the same answer as the real call (what would change) and changes nothing.
- A command writes only through its `tx`, never commits, and does nothing outside the documents; an effect on the world
  is an outbox entry performed after commit (§all-or-nothing). A port without transactions honors `dryRun` itself.
- A mutation that would change nothing returns a coded result or is refused, never a bare success (§no-op-writes-say-so).
- The actor comes from the verifier, never from command input; secrets never travel through inputs, outputs or representations.
- A missing entity uses the entity's namespace (`doc.not_found`), not `command.*`.
- Agent policy: a non-agent actor passes; an agent acting for no one may only `read`; an agent is allowed only if both
  it and its `on` actor are allowed on this target; a `destructive` command needs `confirmed` true.
- Models submit commands as transactions through the replica and never reach around them to a store (§one-path-to-behavior).
- Every catalog entry is namespaced, described, has a valid effect and both schemas.

## Conformance cases

- Catalog: name `doc.rename` → valid; `rename` (no namespace) → invalid; empty description → invalid; effect `write` → invalid.
- `invoke("nope", ...)` → `command.not_found`; one record with code `command.not_found`.
- Register `doc.rename` twice → programmer error.
- Stored `{id: "a", title: "old"}`; `invoke("doc.rename", {id: "a", title: "x"}, dryRun)` → `{id: "a", title: "x"}`;
  stored title is still `old`.
- `invoke("doc.rename", {id: "missing", title: "x"})` → `doc.not_found` with a hint.
- Input that fails the schema → `command.invalid_input`.
- Signal already aborted → nothing runs; record code `cancelled`.
- Successful call → record code `ok`; a command raising an uncoded error → record code `error`.
- Agent with no `on`, `mutate` command → `command.forbidden`; same agent, `read` command its access allows → runs.
- Agent allowed, its `on` actor refused → `command.forbidden`.
- Agent with `on`, both allowed, `destructive` command without `confirmed` → `command.confirmation_required`; with
  `confirmed: true` → runs.
- Person actor under `agentPolicy` → the policy does not refuse.

## Open

- §no-op-writes-say-so requires a coded "unchanged" result, but this pattern declares no code for it, and the example
  `doc.rename` returns plain success when the title is already the current one.
