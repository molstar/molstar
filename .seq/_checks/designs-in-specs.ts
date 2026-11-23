// The designs-in-specs check (the sixth slice's stage 1, `execution/sixth-slice/plan.md`): every agreed design's
// decisions that the code follows are in a package's `_spec`, or say why they stay in the design. A design with a
// decision heading (`## Decision`, `## Decisions`, `## Decided …`) has an `## In the specs` section: a table whose rows
// are `| <decision> | <where> |`, `<where>` either `<spec> › <section>` (the section a heading or a bold lead,
// `- **<section>**`, in that spec, matched by its start) or `design-only: <reason>`. The specs:
//   core      packages/seq/core/_spec/README.md
//   viewer    packages/seq/viewer/_spec/README.md
//   ext-rondo packages/seq/ext-rondo/_spec/README.md
//   repr      packages/core/src/repr/_spec/README.md
// Designs with no decision heading (plans, catalogs, the vision) are not checked. `PENDING` lists the designs not yet
// reconciled when the check was written (2026-11-23); it only shrinks, and stage 1 is done when it is empty.
// `--worklist` prints every design that fails, pending or not.
// Usage: designs-in-specs.ts [root] [--worklist]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, trackedFiles } from "../oz/_checks/lib.ts";

const DESIGNS = ".seq/design/";
const SPECS: Record<string, string> = {
  core: "packages/seq/core/_spec/README.md",
  viewer: "packages/seq/viewer/_spec/README.md",
  "ext-rondo": "packages/seq/ext-rondo/_spec/README.md",
  repr: "packages/core/src/repr/_spec/README.md",
};
/** Designs not reconciled when the check was written (file names without `.md`). Only ever shrinks. */
const PENDING = new Set([
  "across-sections", "addressing", "axis", "changes", "composition", "data-store", "dynamic-nodes", "edge-cases",
  "examples", "examples-and-stories", "extensibility", "extension-openings", "files", "format", "genome-coordinates",
  "genome-formats", "grids-and-rows", "immutable", "language-description", "layout-and-draw-lists", "layouts",
  "locations-and-queries", "maps", "merge", "molecules", "packages", "patches", "ports", "properties",
  "protocol", "region-reads", "rendering", "ring-extension", "security", "seq-data-ports", "sources", "stories",
  "structure-syncing", "testing", "themes", "tooltips", "tree-core",
]);

const argv = process.argv.slice(2);
const WORKLIST = argv.includes("--worklist");
if (WORKLIST) process.argv.splice(process.argv.indexOf("--worklist"), 1);
const { root: ROOT } = args();

function section(text: string, title: RegExp): string | undefined {
  const lines = text.split("\n");
  const from = lines.findIndex((l) => /^## /.test(l) && title.test(l.slice(3)));
  if (from < 0) return undefined;
  let to = lines.findIndex((l, i) => i > from && /^## /.test(l));
  if (to < 0) to = lines.length;
  return lines.slice(from + 1, to).join("\n");
}

const specLeads = new Map<string, string[]>();
function leads(spec: string): string[] {
  let l = specLeads.get(spec);
  if (l) return l;
  const text = readFileSync(join(ROOT, SPECS[spec]), "utf8");
  l = [];
  for (const line of text.split("\n")) {
    const h = /^#+\s+(.*)$/.exec(line);
    if (h) l.push(h[1].trim());
    for (const m of line.matchAll(/\*\*([^*]+)\*\*/g)) l.push(m[1].trim());
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
      if (!m || !(m[1] in SPECS)) {
        own.push(`${f}: "${cells[0]}": where is neither "<spec> › <section>" nor "design-only: <reason>"`);
        continue;
      }
      const want = norm(m[2]);
      if (!leads(m[1]).some((l) => norm(l).startsWith(want)))
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
    problems.push(`${DESIGNS}${id}.md: passes now: drop it from PENDING`);
finish("designs-in-specs", problems);
