// The designs-in-specs check (§mechanical-drift-checks; a scaffold is never cited by code, so what the code follows
// lives in the packages' `_spec`): every agreed design's decisions that the code follows are in a package's `_spec`,
// or say why they stay in the design. A design with a decision heading (`## Decision`, `## Decisions`, `## Decided …`)
// has an `## In the specs` section: a table whose rows are `| <decision> | <where> |`, `<where>` either
// `<spec> › <section>` (the section a heading or a bold lead, `- **<section>**`, in that spec, matched by its start;
// `<spec> › <heading> › <section>` names headings above it, outermost first, each matched by its start) or
// `design-only: <reason>`. Designs with no decision heading (plans, catalogs, the vision) are not checked.
// A pointer must be unambiguous: when its text starts more than one heading or bold lead of the spec, it fails and
// names the candidates with their headings, unless exactly one of them is the pointer's whole text; a heading above
// it in the pointer settles it. Why: the first match was taken, and about ten pointers landed on the wrong section
// while the check passed.
// `--worklist` prints every design that fails, pending or not; `--self-check` resolves planted pointers.
// The method's, not an instance's: the rule holds for every scaffold that governs code. The instance's part is
// `<instance>/_checks/designs-in-specs.json`, required: `{ "specs": { "<key>": "<path of a _spec file>" },
// "pending": ["<design file name without .md>"] }`; `pending` lists the designs not yet reconciled and only shrinks.
// A spec that does not exist yet reads as one with no sections.
// Usage: designs-in-specs.ts [root] --instance .scaffolds/<name> [--worklist] [--self-check]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { finish, instanceFile, methodArgs, section, trackedFiles } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, flags } = methodArgs({ flags: ["--worklist", "--self-check"] });
const WORKLIST = flags.has("--worklist");
const DESIGNS = `${INSTANCE}/design/`;
const { specs: SPECS = {}, pending = [] } = instanceFile<{ specs?: Record<string, string>; pending?: string[] }>(
  ROOT,
  INSTANCE,
  "designs-in-specs.json",
);
/** Designs not reconciled when the check was written. Only ever shrinks. */
const PENDING = new Set(pending);

/** A section of a spec: a heading or a bold lead, with the headings it sits under, outermost first. */
interface Lead {
  text: string;
  under: string[];
}

const specLeads = new Map<string, Lead[]>();
function leads(spec: string): Lead[] {
  let l = specLeads.get(spec);
  if (l) return l;
  const path = join(ROOT, SPECS[spec] as string);
  l = parseLeads(existsSync(path) ? readFileSync(path, "utf8") : "");
  specLeads.set(spec, l);
  return l;
}

/** The headings and bold leads of a spec, each with the headings above it. */
function parseLeads(text: string): Lead[] {
  const out: Lead[] = [];
  const stack: { level: number; text: string }[] = [];
  for (const line of text.split("\n")) {
    const h = /^(#+)\s+(.*)$/.exec(line);
    if (h) {
      const level = (h[1] as string).length;
      while (stack.length && (stack[stack.length - 1]?.level ?? 0) >= level) stack.pop();
      const text = (h[2] as string).trim();
      out.push({ text, under: stack.map((s) => s.text) });
      stack.push({ level, text });
    }
    for (const m of line.replace(/\\\*/g, "∗").matchAll(/\*\*([^*]+)\*\*/g))
      out.push({ text: (m[1] as string).replace(/∗/g, "*").trim(), under: stack.map((s) => s.text) });
  }
  return out;
}
const norm = (s: string) => s.replace(/\\/g, "").replace(/[`*]/g, "").toLowerCase().trim();

/** Whether each of `parents` starts a heading of `under`, in order (not necessarily adjacent). */
function within(under: readonly string[], parents: readonly string[]): boolean {
  let i = 0;
  for (const p of parents) {
    while (i < under.length && !norm(under[i] as string).startsWith(p)) i++;
    if (i++ >= under.length) return false;
  }
  return true;
}

/** Where a pointer (`<section>`, or `<heading> › … › <section>` naming headings above it) lands among a spec's leads:
 * the one section, or why not (none, or the candidates it starts, each with its headings). */
function resolve(specLeads: readonly Lead[], pointer: string): { ok: true } | { ok: false; why: string } {
  const parts = pointer.split("›").map(norm);
  const want = parts.pop() as string;
  const hits = specLeads.filter((l) => norm(l.text).startsWith(want) && within(l.under, parts));
  if (hits.length === 0) return { ok: false, why: "no section" };
  if (hits.length === 1 || hits.filter((l) => norm(l.text) === want).length === 1) return { ok: true };
  const named = hits.map((h) => `"${[...h.under.slice(1), h.text].join(" › ")}"`);
  return { ok: false, why: `${hits.length} sections start with it: ${named.join(", ")}` };
}

if (flags.has("--self-check")) {
  const planted = parseLeads(
    [
      "# Spec",
      "## Language",
      "- **Links**: a",
      "- **Links between rows**: b",
      "- **Bounds of the axis**: c",
      "- **Bounds of the ring**: d",
      "### Configuration",
      "## Resolved view",
      "### Configuration",
      "- **Bridges**: e",
      "## Rendering",
      "### Bridges",
    ].join("\n"),
  );
  const cases: [string, boolean][] = [
    ["Links", true], // the whole text of one, a prefix of another
    ["Links between", true],
    ["Bounds", false], // starts two, the whole text of neither
    ["Bounds of the ring", true],
    ["Configuration", false], // two sections of one name
    ["Language › Configuration", true], // settled by the heading above it
    ["Bridges", false],
    ["Resolved view › Bridges", true],
    ["Rendering › Bridges", true],
    ["Rendering › Links", false], // not under that heading
    ["Ranges", false], // none
  ];
  const wrong = cases
    .filter(([p, ok]) => resolve(planted, p).ok !== ok)
    .map(([p, ok]) => `"${p}" resolved ${ok ? "ambiguous or missing" : "to one section"}, expected the opposite`);
  if (!wrong.length) console.log(`designs-in-specs --self-check: ${cases.length} planted pointers resolved as expected`);
  finish("designs-in-specs --self-check", wrong);
}

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
      const at = resolve(leads(m[1] as string), m[2] as string);
      if (!at.ok)
        own.push(
          at.why === "no section"
            ? `${f}: "${cells[0]}": ${m[1]} has no section "${m[2]}"`
            : `${f}: "${cells[0]}": "${m[2]}" is ambiguous in ${m[1]}, ${at.why}: write more of its text, or a heading above it (<spec> › <heading> › <section>)`,
        );
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
