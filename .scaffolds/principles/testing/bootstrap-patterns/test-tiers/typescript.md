# Test tiers and targeted runs: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Naming = tier

```
src/foo/foo.quick.test.ts          # unit, in-memory, ms
src/foo/foo.medium.test.ts         # module compositions, in-memory
src/foo/foo.guard.test.ts          # allocation / heap / scaling budgets (§performance)
src/foo/foo.full.test.ts           # slow or external (network, real provider)
src/foo/foo.eval.ts                # a real model or other metered service: only with acknowledgement and a budget (§evals)
{apps,services}/<product>/<name>/smoke/*.smoke.test.ts   # the built artifact, every target it ships, no external services
```

## Runner (vitest example)

```ts
// tooling/scripts/test.ts
// usage: tsx tooling/scripts/test.ts <quick|medium|guard|smoke|full|gate|nightly> [dir]
import { spawnSync } from "node:child_process";

const TIERS = { quick: ".quick.test", medium: ".medium.test", guard: ".guard.test", smoke: ".smoke.test", full: ".full.test" } as const;
type Tier = keyof typeof TIERS;
const SETS: Record<string, readonly Tier[]> = {
  gate:    ["quick", "medium", "guard", "smoke"],              // unfiltered, before every merge (§the-gate)
  nightly: ["quick", "medium", "guard", "smoke", "full"],
};

const [what = "", dir] = process.argv.slice(2);
const tiers = Object.hasOwn(SETS, what) ? SETS[what]! : Object.hasOwn(TIERS, what) ? [what as Tier] : [];
if (!tiers.length) { console.error(`unknown tier or set: "${what}"`); process.exit(2); }
// A directory makes it a *targeted* run: for development, never for the gate or nightly.
// Evals are not a tier of this runner: `mise run evals -- <path> --ack --runs N --budget <usd>` refuses to start without all three,
// stops when the budget is spent, and writes results under the gitignored .tmp/ (§evals, §machine-invariant).
if (dir && Object.hasOwn(SETS, what)) { console.error("gate and nightly run unfiltered"); process.exit(2); }

// vitest's positional filters match path substrings and are OR'd: scope by --dir, select tiers by suffix.
const args = ["vitest", "run", ...(dir ? ["--dir", dir] : []), ...tiers.map((t) => TIERS[t])];
process.exit(spawnSync("npx", args, { stdio: "inherit" }).status ?? 1);
```

Typical agent loop:

```bash
mise run test -- quick lib/storage                 # targeted, first
mise run test -- quick                             # then widen
mise run gate                                      # before merge: build + unfiltered tiers + checks
```

Tiers, sets and the time rule: see [README.md](README.md#invariants).

## Conformance suites: the file's tier bounds the cases

Shared behavioral checks use the conformance library ([§conformance](../conformance/README.md)). Cases carry a cost tag, and a registering file runs cases only up to its own tier.

```ts
// memory-store.quick.test.ts
register(documentStoreSuite, () => memoryStore(), { label: "memory", upTo: "quick" });
// cloud-store.full.test.ts: needs a real deployment, so full tier, all cases
register(documentStoreSuite, () => cloudStore(testClient(), uniqueTable()), { label: "cloud", upTo: "full" });
```

## Composition and smoke tests for the final app

```ts
// packages/product-x/compose/compose.medium.test.ts: real wiring (§ports-and-adapters), in-memory adapters, driven through commands (§commands)
test("create, rename, observe", async () => {
  const scope = new Scope();
  const app = createApp(local(scope), scope);
  const { id } = (await app.invoke("doc.create", { title: "a" })) as { id: string };
  await app.invoke("doc.rename", { id, title: "b" });
  expect(await app.invoke("doc.get", { id })).toMatchObject({ title: "b" });   // observe through a read command (§reachable-queries)
  await scope.dispose();
});
```

```ts
// services/product-x/api/smoke/boot.smoke.test.ts: the BUILT artifact on the local composition, no external services (§whole-composition-tested)
import { spawn, type ChildProcess } from "node:child_process";
import { fileURLToPath } from "node:url";
declare function readPort(proc: ChildProcess): Promise<number>;   // resolves on the "listening on <port>" line

const invoke = (port: number, name: string, input: unknown) =>
  fetch(`http://localhost:${port}/commands/${name}`, {
    method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(input),
  });

test("built app boots and answers one read command per area", async () => {
  const server = fileURLToPath(new URL("../dist/server.js", import.meta.url));   // vitest's cwd is the repo root
  const proc = spawn("node", [server], { env: { ...process.env, COMPOSITION: "local", PORT: "0" } });
  try {
    const port = await readPort(proc);
    for (const name of ["doc.list" /* , one read command per major area */]) {   // the project's read commands (§commands)
      expect((await invoke(port, name, {})).ok, name).toBe(true);
    }
  } finally { proc.kill(); }
}, 30_000);
```

Smokes, provider-backed checks and waiting in tests: see [README.md](README.md#shape). In vitest, wait with
`await vi.waitFor(() => expect(x).toBe(y))`.

## Guard tier: allocation and complexity

Guard tests run node with `--expose-gc` (for vitest: `poolOptions.forks.execArgv`).

```ts
// render.guard.test.ts
const time = (f: () => void) => { const t = performance.now(); f(); return performance.now() - t; };

it("steady-state frame allocates nothing", () => {
  if (!global.gc) throw new Error("guard tier needs --expose-gc");
  const r = createRenderer(preallocated());
  for (let i = 0; i < 100; i++) r.frame();       // warm up: JIT, lazy init
  global.gc();                                   // collect BEFORE measuring...
  const before = process.memoryUsage().heapUsed;
  for (let i = 0; i < 1_000; i++) r.frame();     // ...and not after: garbage made here must still be on the heap.
  expect(process.memoryUsage().heapUsed - before).toBeLessThan(16 * 1024);   // keep the loop well under the young generation
});

it("index build scales ~linearly", () => {
  const best = (n: number) => {
    const items = makeItems(n);
    let t = Infinity;
    for (let i = 0; i < 5; i++) t = Math.min(t, time(() => buildIndex(items)));   // early runs warm the JIT; keep the best
    return t;
  };
  expect(best(40_000) / best(10_000)).toBeLessThan(8);   // 4x input; quadratic would be ~16
});
```

Helpers assumed: `createRenderer`, `preallocated`, `makeItems`, `buildIndex`.
