# Spec mirror {#spec-mirror}

Rules §spec-folders, §mechanical-drift-checks: `_spec/` is a short semantic mirror of the code as it is, linking to the registries instead of copying them; spec = map, code = territory; one gate check keeps links resolving and flags code that changed without its entry.

## Layout

```
_spec/                               # global: architecture, glossary, ports index
  architecture.md
  ports.md                           # every port -> adapters -> where composed
  glossary.md
lib/storage/
  AGENTS.md                          # how to work here: scripts, local conventions
  _spec/
    module.md                        # purpose, public surface, invariants, non-goals
    memory-store.md                  # one file per significant unit
```

## Spec entry shape

```md
---
code: [lib/storage/src/memory.ts]   # path(s) this entry describes
kind: adapter
implements: DocumentStore
---
## Purpose
In-memory DocumentStore for tests, stories and local evals.

## Invariants
- `watch(id)` delivers the current value first, then every change.
- `get` of a missing id resolves to `undefined`; it never throws.

## Links
- Port: [module.md](module.md#documentstore)
- Conformance: `suite:documentStoreSuite`
- Errors: `error:store.conflict`
```

The entry states meaning (purpose, invariants, links); signatures, command lists and error catalogs are linked by registry reference, never copied.

## The gate check

Part of the drift check run at the gate (§mechanical-drift-checks; CI, see the testing module, if used). Broken links and unknown registry references block. "Code changed but its entry did not" is advice: often the entry is still true.

```ts
// tooling/scripts/drift-check.ts (spec half; the drift-note half is in finalize, §drift-triage)
import { execFileSync } from "node:child_process";
import { existsSync, readFileSync, readdirSync } from "node:fs";
import { dirname, join } from "node:path";
import matter from "gray-matter";
declare function registryNames(): Iterable<string>;   // "command:doc.rename", "error:store.conflict", "suite:...", built from the registries

const git = (...args: string[]) => execFileSync("git", args, { encoding: "utf8" }).trim();
const base = git("merge-base", "HEAD", process.env.BASE_REF ?? "origin/main");
const changed = git("diff", "--name-only", base, "HEAD").split("\n").filter(Boolean);   // changed on this branch

function* specFiles(dir: string): Generator<string> {
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    if (e.name === "node_modules" || (e.name.startsWith(".") && e.name !== "_spec")) continue;   // skip .git, .execution, ...
    const p = join(dir, e.name);
    if (e.isDirectory()) yield* specFiles(p);
    else if (p.endsWith(".md") && p.split(/[\\/]/).includes("_spec")) yield p;
  }
}

const known = new Set(registryNames());
const errors: string[] = [];
for (const f of specFiles(".")) {
  const { data, content } = matter(readFileSync(f, "utf8"));
  const code = ([] as string[]).concat(data.code ?? []);
  if (!changed.includes(f) && code.some((c) => changed.some((p) => p === c || p.startsWith(`${c}/`))))
    console.warn(`advice: ${f}: its code changed on this branch; check the entry`);
  for (const [, href] of content.matchAll(/\]\(([^)#\s]+)[^)]*\)/g))
    if (!/^[a-z]+:/.test(href) && !existsSync(join(dirname(f), href))) errors.push(`${f}: broken link ${href}`);
  for (const [, ref] of content.matchAll(/`((?:command|error|suite|concept):[\w.-]+)`/g))
    if (!known.has(ref)) errors.push(`${f}: unknown registry reference ${ref}`);
}
errors.forEach((e) => console.error(e));
process.exit(errors.length ? 1 : 0);
```

An agent that changes code updates the matching entry in the same change, or files a followup issue if the gap cannot be closed there (§mechanical-drift-checks).

## AGENTS.md skeleton

```md
# storage
Principles: .principles/README.md; here mainly the architecture module (structure, ports, API boundaries, composition) and the testing module, if used.
Tests: the entry point's `test` task, quick tier, `lib/storage` (`mise run test -- quick lib/storage` in the examples)
Ports defined here: DocumentStore (adapters: memory here; a provider's in lib/storage-<provider>-runtime)
Local rules: no provider imports; every adapter passes documentStoreSuite.
Spec: ./_spec/module.md
```
