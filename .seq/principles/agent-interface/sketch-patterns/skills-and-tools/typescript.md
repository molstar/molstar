# Skills and tools: the agent interface: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## The whole tool surface

```ts
// generated once at startup from the list of domains; it never changes shape afterwards (§fixed-tool-shape)
export const DOMAINS = ["doc", "scene", "mol", "data", "sem", "agent", "wf"] as const;   // agent: 23, wf: 24

const exprInput = {
  type: "object",
  properties: { expr: { type: "string", description: "One form (S-expression, or the JSON spelling of the same tree). Several forms: (do f1 f2 ...) to read, (tx ...) to write." } },
  required: ["expr"],
} as const;

export const tools = [
  { name: "skills_read", description: "Skill index, content and search: (index), (skill <domain> [section]), (search \"text\" [k]).", inputSchema: exprInput },
  ...DOMAINS.flatMap((d) => [
    { name: `${d}_read`,  description: `Read-only ${d} access: queries, (get <locator>), (help ...), (more <cursor>). Never mutates.`, inputSchema: exprInput },
    { name: `${d}_write`, description: `Change ${d}: commands, (dry-run <form>), (tx ...) as one transaction. Subject to policy.`, inputSchema: exprInput },
  ]),
];
```

Why the surface is fixed: see [README.md](README.md#invariants).

## Sessions

```ts
/** Per agent session: identity, policy, cursors for paged results. Created by the transport adapter (MCP, CLI, tests). */
export interface ToolSession {
  actor: { kind: "agent"; id: string };
  policy: Policy;                                        // 08: least privilege, confirmation (§least-privilege)
  cursors: Map<string, { items: readonly Node[]; next: number }>;
  signal?: AbortSignal;                                  // the transport's cancellation
  budget: { items: number; chars: number };              // result budget per call (§compact-results)
  /** Set by the transport after the USER approved this exact form (§least-privilege); never read from the agent's own input. */
  approved?: (form: Node) => boolean;
}
```

Where approvals come from: see [README.md](README.md#shape).

## Read dispatch

```ts
export function domainRead(domain: string, deps: { langs: LanguageRegistry; registry: CommandRegistry; state: TreeState; ctx: QueryContext }) {
  return async (expr: string, s: ToolSession): Promise<string> => {
    const { node, repairs } = readExpr(expr);                                  // lenient + JSON spelling (below)
    const root = deps.state.value;                                             // one snapshot for the whole call, incl. (do ...)
    const ctx = { ...deps.ctx, signal: s.signal };

    const one = async (f: Node): Promise<Node> => {
      switch (tag(f)) {
        case "help": return help(domain, deps, (f as readonly Node[]).slice(1));          // dynamic skill content (§dynamic-skills-follow-system)
        case "get":  return getAt(root, toLocator((f as readonly Node[])[1]!));             // point read by locator (§semantic-layer)
        case "more": return more(s, nameOf((f as readonly Node[])[1]));
        case "do":   return [sym("results"), ...(await Promise.all((f as readonly Node[]).slice(1).map(one)))];
      }
      const cmd = deps.registry.find(tag(f) ?? "");
      if (cmd) {                                                               // read-effect commands only
        if (cmd.effect !== "read") throw toolErrors.error("tool.wrong_tool", undefined, { hint: `use ${domain}_write` });
        return toNode(await deps.registry.invoke(cmd.name, formToArgs(f), { dryRun: false, actor: s.actor, signal: s.signal }, s.policy));
      }
      const { lang, canonical } = await resolveQuery(deps.langs, domain, print(f));          // any of the domain's query languages (§queryable-representations)
      const paths = await deps.langs.engine(domain, lang.id).runAsync(root, canonical, ctx);
      return page(s, paths.map((p) => locatorAt(domain, root, p)));                          // compact, with a cursor if large
    };

    return print(withRepairs(await one(node), repairs));
  };
}
```

Command scoping: see [Open](README.md#open).

## Write dispatch

```ts
export function domainWrite(domain: string, deps: { registry: CommandRegistry; runTx: (tx: Tx, ctx: TxCtx) => Promise<TxResult> }) {
  return async (expr: string, s: ToolSession): Promise<string> => {
    const { node, repairs } = readExpr(expr);
    const run = (f: Node, dryRun: boolean) => {
      const name = tag(f) ?? "";
      if (!deps.registry.find(name)) throw commandError("command.not_found", `no command ${name}`, { hint: `(help commands) in ${domain}_read` });
      return deps.registry.invoke(name, formToArgs(f), { dryRun, actor: s.actor, signal: s.signal, confirmed: s.approved?.(f) === true }, s.policy);
    };

    if (tag(node) === "dry-run") return print(withRepairs(toNode(await run((node as readonly Node[])[1]!, true)), repairs));
    if (tag(node) === "tx") {                                                // all or nothing, one commit per document (§batched-calls, §transactions)
      const tx = parseTx(node, { origin: "agent" });                         // (step [name] <form>) ..., (use name field)
      const { results } = await deps.runTx(tx, { actor: s.actor, signal: s.signal });
      return print(withRepairs([sym("results"), ...Object.values(results).map(toNode)], repairs));
    }
    return print(withRepairs(toNode(await run(node, false)), repairs));
  };
}
// Commands write into the transaction's scope (CommandContext.tx, §transactions) and never commit themselves, so (tx ...) is atomic.
// Approval of a destructive step inside a tx: the whole tx yields once with every step that needs it (agent-harness's `harness-yields` (not adopted here)).
```

Refusals and outline records: see [README.md](README.md#shape). The record shape:

```ts
/** One outline record; `can` is policy-filtered for the requester (§one-policy-port), arguments prefilled, holes as `_`. */
export interface OutlineRecord { readonly kind: string; readonly at: Locator; readonly excerpt: string;
  readonly status?: string; readonly children?: number; readonly hash: string; readonly can: readonly Node[] }
// `can` is computed from the concept's commands (meso's `documents` (not adopted here)) and `policy.allows(requester, cmd, at)`; an empty list is an answer
```

## Forms, JSON spelling, results

```ts
/** One input string -> one tree. S-expression by default; JSON is another spelling of the same tree (§one-recursive-form). */
export function readExpr(text: string): { node: Node; repairs: readonly Repair[] } {
  const t = text.trim();
  if (t.startsWith("{") || t.startsWith("[")) return { node: fromShorthand(JSON.parse(t)), repairs: [] };
  const { node, repairs } = reprArg(t);                                      // lenient, reported (§commands, §parser-repairs-slips)
  return { node, repairs };
}

/** {"doc.rename": {"id": "a", "title": "b"}}  ->  (doc.rename (id "a") (title "b")); arrays are lists.
 *  The agent convenience spelling (meso's `convenience-spellings` (not adopted here)), NOT 21's lossless `fromJson` ($sym/$big tags). Open: accept both at the edge? */
export const fromShorthand = (v: unknown): Node =>
  Array.isArray(v) ? v.map(fromShorthand)
  : v && typeof v === "object"
    ? (([k, body]) => [sym(k), ...Object.entries(body as object).map(([a, x]) => [sym(a), fromShorthand(x)])])(Object.entries(v)[0]!)
  : (v as Atom);

/** (doc.rename (id "a") (title "b")) -> { id: "a", title: "b" }; an attribute with several values becomes an array
 *  ((steps (step ..) (step ..))). Open: single vs list for one value is left to the command's schema (§commands). */
export const formToArgs = (f: Node): Record<string, unknown> =>
  Object.fromEntries((f as readonly Node[]).slice(1).filter(Array.isArray).map((a) => { const xs = a as readonly Node[];
    return [nameOf(xs[0]), xs.length > 2 ? xs.slice(1).map(toPlain) : toPlain(xs[1])]; }));

/** Results over budget become a page with a cursor, never a dump (§compact-results). */
export function page(s: ToolSession, items: readonly Node[]): Node {
  if (items.length <= s.budget.items) return [sym("items"), ...items];
  const id = newId("c");
  s.cursors.set(id, { items, next: s.budget.items });
  return [sym("page"), [sym("count"), items.length], [sym("items"), ...items.slice(0, s.budget.items)], [sym("more"), id]];
}

export const toolErrors = defineErrors("tool", {
  "tool.wrong_tool":      { description: "This form mutates; send it to <domain>_write.", retryable: false },
  "tool.unknown_cursor":  { description: "Cursor expired or unknown; rerun the query.", retryable: false },
});
// helpers assumed: getAt, more (throws tool.unknown_cursor), locatorAt and toLocator (§semantic-layer), toNode, toPlain,
// withRepairs (adds (normalized ...) (repairs ...) when repairs exist)
```

A session, end to end: see [README.md](README.md#shape).

## Skills: static files plus generated content

The skill files and a sample `SKILL.md`: see [README.md](README.md#shape). The generated part:

```ts
/** (help ...) answers from the registries of the RUNNING system, so skills cannot drift (§dynamic-skills-follow-system). */
export function help(domain: string, deps: { langs: LanguageRegistry; registry: CommandRegistry }, args: readonly Node[]): Node {
  switch (args.length ? nameOf(args[0]) : "topics") {
    case "commands":  return [sym("commands"), ...deps.registry.catalog().filter((c) => c.name.startsWith(`${domain}.`)).map((c) => [sym(c.name), c.effect, c.description])];
    case "functions": return [sym("functions"), ...(deps.langs.get(domain).queries ?? []).flatMap((q) =>
                        functionsOf(q.lang).map((f) => [sym(f.name), sym(f.role), f.doc]))];      // every query language of the domain
    case "dialects":  return [sym("dialects"), ...(deps.langs.get(domain).dialects ?? []).map((d) => [sym(d.id), d.lowersTo])];
    case "errors":    return [sym("errors"), ...errorCatalogFor(domain)];
    default:          return [sym("topics"), sym("commands"), sym("functions"), sym("dialects"), sym("errors"), sym("examples")];
  }
}
// helpers assumed: functionsOf (below), errorCatalogFor (the domain's defineErrors specs as nodes, §error-codes)
```

## Search: fuzzy now, semantic later

```ts
/** One entry per thing an agent might look for, derived from skills and registries (never hand-maintained). */
export interface SearchDoc {
  loc: Node;                    // (skill mol "Workflows"), (command mol.set-color), (function mol.select atoms), (code repr.link.cycle)
  title: string;                // the name an agent would type
  text: string;                 // description, synonyms, doc string
}

/** The search port: lexical today; an embedding backend implements the same interface later (§searchable-skills). */
export interface SkillSearch {
  search(text: string, k: number, signal?: AbortSignal): Promise<readonly { loc: Node; title: string; score: number }[]>;
}

export function buildSearchDocs(skills: SkillFiles, registry: CommandRegistry, langs: LanguageRegistry): SearchDoc[] {
  return [
    ...skills.sections().map((s) => ({ loc: [sym("skill"), s.domain, s.heading], title: s.heading, text: s.body })),
    ...registry.catalog().map((c) => ({ loc: [sym("command"), sym(c.name)], title: c.name, text: c.description })),
    ...langs.all().flatMap((l) => (l.queries ?? []).flatMap((q) => functionsOf(q.lang).map((f) =>
      ({ loc: [sym("function"), sym(q.lang.id), sym(f.name)], title: f.name, text: `${f.role}: ${f.doc}` })))),
    ...allErrorCatalogs().map((e) => ({ loc: [sym("code"), sym(e.code)], title: e.code, text: e.description })),
  ];
}

/** Default: fuzzy lexical. Titles weigh more than text; trigram similarity tolerates typos and partial names. */
export function fuzzySearch(docs: readonly SearchDoc[]): SkillSearch {
  const index = docs.map((d) => ({ d, t: trigrams(d.title), x: trigrams(d.text) }));    // built once per registry version
  return {
    async search(text, k) {
      const q = trigrams(text);
      return index
        .map(({ d, t, x }) => ({ loc: d.loc, title: d.title, score: 2 * jaccard(q, t) + jaccard(q, x) }))
        .filter((r) => r.score > 0.1)
        .sort((a, b) => b.score - a.score)
        .slice(0, k);                                          // N is hundreds to low thousands: a scan is fine (§efficient-structures)
    },
  };
}

/** Later: same port, embeddings computed when the registries change; can be combined with fuzzy (hybrid ranking). */
export declare function semanticSearch(docs: readonly SearchDoc[], embed: EmbeddingPort): SkillSearch;
// helpers assumed: functionsOf (name, role, doc per registered function, §query-dsls), trigrams, jaccard, allErrorCatalogs, EmbeddingPort (a runtime port bound in the composition root)
```

Search example and what `(skill <domain>)` returns: see [README.md](README.md#shape). The composition root binds `fuzzySearch` (or later `semanticSearch`); `skills_read` and the tool list stay the same.

## Conformance

```ts
export const toolSurfaceSuite = suite<ToolHost>("ToolSurface",
  kase("tool list equals the golden list (names and input schema never change)", "quick", (h) => { /* snapshot */ }),
  kase("every catalog command is reachable through <domain>_read or _write", "quick", (h) => { /* ... */ }),
  kase("_read never mutates: document hash unchanged after any read form", "quick", async (h) => { /* ... */ }),
  kase("results stay within budget; large results page with a cursor", "quick", async (h) => { /* ... */ }),
  kase("JSON spelling and S-expression spelling give the same result", "quick", async (h) => { /* ... */ }),
  kase("(tx ...) on _write is one commit per document, or none when any step fails", "medium", async (h) => { /* ... */ }),
);

export const skillSuite = suite<SkillFiles>("Skills",
  kase("every example form runs through the real tools and matches its expected output", "medium", async (k) => { /* ... */ }),
  kase("every command, function and dialect named in a skill exists", "quick", (k) => { /* ... */ }),
  kase("each skill fits its context budget", "quick", (k) => { /* e.g. under 1500 tokens */ }),
  kase("search finds every command, function and code by its exact name, and by a one-typo variant", "quick", (k) => { /* ... */ }),
);
```
