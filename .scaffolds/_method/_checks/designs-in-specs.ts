// The designs-in-specs check (§mechanical-drift-checks; a scaffold is never cited by code, so what the code follows
// lives in the packages' `_spec`): every agreed design's decisions that the code follows are in a package's `_spec`,
// or say why they stay in the design. A design with a decision heading (`## Decision`, `## Decisions`, `## Decided …`)
// has an `## In the specs` section: a table whose rows are `| <decision> | <where> |`, `<where>` either
// `<spec> › <section>` (the section a heading or a bold lead, `- **<section>**`, in that spec, matched by its start) or
// `design-only: <reason>`. Designs with no decision heading (plans, catalogs, the vision) are not checked.
// `--worklist` prints every design that fails, pending or not.
// The method's, not an instance's: the rule holds for every scaffold that governs code. The instance's part is
// `<instance>/_checks/designs-in-specs.json`, required: `{ "specs": { "<key>": "<path of a _spec file>" },
// "pending": ["<design file name without .md>"] }`; `pending` lists the designs not yet reconciled and only shrinks.
// A spec that does not exist yet reads as one with no sections.
// Usage: designs-in-specs.ts [root] --instance .scaffolds/<name> [--worklist]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { finish, instanceFile, methodArgs, section, trackedFiles } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, flags } = methodArgs({ flags: ["--worklist"] });
const WORKLIST = flags.has("--worklist");
const DESIGNS = `${INSTANCE}/design/`;
const { specs: SPECS = {}, pending = [] } = instanceFile<{ specs?: Record<string, string>; pending?: string[] }>(
  ROOT,
  INSTANCE,
  "designs-in-specs.json",
);
/** Designs not reconciled when the check was written. Only ever shrinks. */
const PENDING = new Set(pending);

const specLeads = new Map<string, string[]>();
function leads(spec: string): string[] {
  let l = specLeads.get(spec);
  if (l) return l;
  const path = join(ROOT, SPECS[spec] as string);
  const text = existsSync(path) ? readFileSync(path, "utf8") : "";
  l = [];
  for (const line of text.split("\n")) {
    const h = /^#+\s+(.*)$/.exec(line);
    if (h) l.push((h[1] as string).trim());
    for (const m of line.replace(/\\\*/g, "∗").matchAll(/\*\*([^*]+)\*\*/g))
      l.push((m[1] as string).replace(/∗/g, "*").trim());
  }
  specLeads.set(spec, l);
  return l;
}
const norm = (s: string) => s.replace(/\\/g, "").replace(/[`*]/g, "").toLowerCase().trim();

const problems: string[] = [];
const worklist: string[] = [];
for (const f of trackedFiles(ROOT)) {
  if (!f.startsWith(DESIGNS) || !f.endsWith(".md") || !existsSync(join(ROOT, f))) continue;
  const id = f.slice(DESIGNS.length, -3);
  const text = readFileSync(join(ROOT, f), "utf8");
  if (section(text, /^Decisions?\b|^Decided\b/) === undefined) continue;
  const own: string[] = [];
  const table = section(text, /^In the specs\b/);
  if (table === undefined) own.push(`${f}: no "## In the specs" section`);
  else {
    const rows = table.split("\n").filter((l) => /^\|/.test(l) && !/^\|\s*-/.test(l)).slice(1);
    if (rows.length === 0) own.push(`${f}: "## In the specs" has no rows`);
    for (const row of rows) {
      const cells = row.split("|").slice(1, -1).map((c) => c.trim());
      const where = cells[1] ?? "";
      if (/^design-only:\s*\S/.test(where)) continue;
      const m = /^`?([a-z-]+)`?\s*›\s*(.+)$/.exec(where);
      if (!m || !((m[1] as string) in SPECS)) {
        own.push(`${f}: "${cells[0]}": where is neither "<spec> › <section>" nor "design-only: <reason>"`);
        continue;
      }
      const want = norm(m[2] as string);
      if (!leads(m[1] as string).some((l) => norm(l).startsWith(want)))
        own.push(`${f}: "${cells[0]}": ${m[1]} has no section "${m[2]}"`);
    }
  }
  if (own.length === 0) continue;
  worklist.push(...own);
  if (!PENDING.has(id)) problems.push(...own);
}
if (WORKLIST) console.log(`designs-in-specs: owed\n${worklist.map((w) => `  - ${w}`).join("\n")}`);
for (const id of PENDING)
  if (!worklist.some((w) => w.startsWith(`${DESIGNS}${id}.md:`)) && existsSync(join(ROOT, `${DESIGNS}${id}.md`)))
    problems.push(`${DESIGNS}${id}.md: passes now: drop it from designs-in-specs.json's pending`);
finish("designs-in-specs", problems);
