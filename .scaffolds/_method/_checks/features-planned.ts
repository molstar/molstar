// The features-planned check (the rule "the features table moves with the merges", from Seq*'s fold-backs): a row of
// the scaffold's features table whose status says "planned" or "to come" fails once every lane of its design is merged
// in the plan (the lane that builds a feature moves its row in its own change). The rows name designs, not lanes: the
// plan's "The designs" table (`| [<design>](<file>.md) ... | <lane>, <lane>; ... |`, the method's plan template) maps
// each design to its lanes, and a lane is merged by a line `- [x] **<LANE>**`, `- [x] **<LANE>, <LANE>**` or
// `- [x] <LANE>` in the `checklist.md` beside the plan (§plan-and-checklist), or in the plan itself (the folders
// written before the checklist). A row names its designs by the links of its Design column; a row whose designs are
// not in the table cannot be checked. The status column is found by the header's index (`Feature`, `Status`, `Design`).
// The method's, not an instance's: it reads only the method's plan shape and a features table; a scaffold with no
// features table passes `--features` to another or does not run it. No instance file: the plan is always given, since
// a check never names a work folder (§execution-unreferenced).
// The same file holds the brief check (`--brief`, "items landed, not reported", loop.md Each step 5): a lane's brief
// numbers its items (`1. **<name>**: ...`), and the lane's row in `findings.md` (by the brief's lane id, `# Lane <ID>:`)
// answers each in its Items cell, `1 ✓ <proof>; 2 ✗ not landed: <where it went>`. An item neither ticked with a proof
// nor marked "not landed", or a lane with no row, fails; each is named. It reads only the method's brief and findings
// shapes. A correction is a new row (findings are append-only), so the lane's last row counts.
// Usage: features-planned.ts [root] --instance .scaffolds/<name> --plan <plan.md, from the root>
//        [--features <table, from the root; default <instance>/design/features.md>]
//        | --brief <brief.md> --findings <findings.md> (from the root) | [--self-check]
//        (self-check: planted plans, checklists, tables, briefs and rows, each failing or passing as it should)
import { existsSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { finish, methodArgs, usage } from "./lib.ts";

const {
  root: ROOT,
  instance: INSTANCE,
  values,
  flags,
} = methodArgs({ values: ["--plan", "--features", "--brief", "--findings"], flags: ["--self-check"] });

const LaneId = /\b[A-Z][A-Z0-9]+\b/g;

/** The planned rows whose designs' lanes are all merged: by the plan's table, and the ticks of the checklist beside
 * it, else (a folder written before checklists) the plan's own. */
function judge(planText: string, checklistText: string | undefined, features: string, where: string): string[] {
  const lanesOf = new Map<string, string[]>();
  for (const line of planText.split("\n")) {
    if (!line.startsWith("|")) continue;
    const cells = line.split("|").slice(1, -1).map((c) => c.trim());
    const design = /^\[([^\]]+)\]\(/.exec(cells[0] ?? "")?.[1];
    if (design && cells[1]) lanesOf.set(design, cells[1].match(LaneId) ?? []);
  }
  const merged = new Set<string>();
  for (const line of `${planText}\n${checklistText ?? ""}`.split("\n")) {
    const m = /^\s*- \[x\]\s+\**((?:[A-Z][A-Z0-9]*)(?:\s*[,;]\s*[A-Z][A-Z0-9]*)*)/.exec(line);
    if (m) for (const id of (m[1] as string).match(LaneId) ?? []) merged.add(id);
  }
  const problems: string[] = [];
  let header: string[] = [];
  for (const line of features.split("\n")) {
    if (!line.startsWith("|") || line.startsWith("|---")) continue;
    const cells = line.split("|").slice(1, -1).map((c) => c.trim());
    if (cells[0] === "Feature") {
      header = cells;
      continue;
    }
    const status = cells[header.indexOf("Status")] ?? "";
    if (!/\b(planned|to come)\b/.test(status)) continue;
    const designs = [...(cells[header.indexOf("Design")] ?? "").matchAll(/\[[^\]]*\]\(([^)]+?)(?:\.md)?(?:#[^)]*)?\)/g)]
      .map((m) => m[1] as string)
      .filter((d) => lanesOf.has(d));
    const lanes = designs.flatMap((d) => lanesOf.get(d) ?? []);
    if (lanes.length && lanes.every((l) => merged.has(l)))
      problems.push(
        `${cells[0]}: status "${status.slice(0, 40)}" but every lane of ${designs.join(", ")} (${lanes.join(", ")}) is merged in ${where} or its checklist.md`,
      );
  }
  return problems;
}

/** The brief's numbered items (`1. **<name>**`, between `**Items**` and the next bold heading line), by number. */
function briefItems(brief: string): Map<number, string> {
  const items = new Map<number, string>();
  let inItems = false;
  for (const line of brief.split("\n")) {
    if (/^\*\*Items\*\*/.test(line)) inItems = true;
    else if (/^\*\*[A-Za-z][^*]*\*\*/.test(line)) inItems = false;
    const m = inItems ? /^(\d+)\.\s+\*\*(.+?)\*\*/.exec(line) : null;
    if (m) items.set(Number(m[1]), m[2] as string);
  }
  return items;
}

/** The items of a brief that a lane's findings row leaves unanswered: each needs `<n> ✓ <proof>` (a proof, not just a
 * tick) or `<n> ✗ not landed...` in the Items cell of the lane's last row. A missing row or no items is a problem. */
function judgeBrief(brief: string, findings: string): string[] {
  const lane = /^# Lane\s+([A-Z][A-Z0-9]*)\b/m.exec(brief)?.[1];
  if (!lane) return ["the brief has no `# Lane <ID>` heading"];
  const items = briefItems(brief);
  if (!items.size) return [`${lane}: the brief has no numbered items (\`1. **<name>**\`) under **Items**`];
  let header: string[] = [];
  let itemsCell: string | undefined;
  for (const line of findings.split("\n")) {
    if (!line.startsWith("|") || line.startsWith("|---")) continue;
    const cells = line.split("|").slice(1, -1).map((c) => c.trim());
    if (cells[0] === "Lane or gate") {
      header = cells;
      continue;
    }
    if (cells[0] === lane) itemsCell = cells[header.indexOf("Items")] ?? "";
  }
  if (itemsCell === undefined) return [`${lane}: no row in findings.md (its first cell is the lane id)`];
  const answers = new Map<number, { mark: string; text: string }>();
  // answers split where the next one starts (`; <n> ✓|✗`), so a proof may hold a semicolon of its own
  for (const part of itemsCell.split(/;\s*(?=\d+\s*[✓✗])/)) {
    const m = /^(\d+)\s*([✓✗])\s*(.*)$/s.exec(part.trim());
    if (m) answers.set(Number(m[1]), { mark: m[2] as string, text: (m[3] as string).trim() });
  }
  const problems: string[] = [];
  for (const [n, name] of items) {
    const a = answers.get(n);
    if (!a) problems.push(`${lane} item ${n} (${name}): neither ticked with a proof nor marked "not landed" in its row`);
    else if (a.mark === "✓" && !a.text) problems.push(`${lane} item ${n} (${name}): ticked without a proof`);
    else if (a.mark === "✗" && !/not landed/i.test(a.text))
      problems.push(`${lane} item ${n} (${name}): marked ✗ without "not landed: <where it went>"`);
  }
  return problems;
}

if (flags.has("--self-check")) {
  const plan = "| Design | Lanes |\n|---|---|\n| [d](../../design/d.md) | AB, CD |\n";
  const table = (status: string) => `| Feature | Status | Design |\n|---|---|---|\n| F | ${status} | [d](d.md) |\n`;
  const cases: [string, string, string | undefined, string, number][] = [
    ["the checklist ticks every lane", plan, "- [x] **AB** (abc)\n- [x] **CD** (def)\n", table("planned"), 1],
    ["the checklist leaves one lane unticked", plan, "- [x] **AB** (abc)\n- [ ] **CD**\n", table("planned"), 0],
    ["an older folder: the plan's own ticks", `${plan}\n- [x] **AB, CD** (abc)\n`, undefined, table("to come"), 1],
    ["an older folder: one lane unticked", `${plan}\n- [x] **AB**\n- [ ] **CD**\n`, undefined, table("planned"), 0],
    ["the row moved with the merge", plan, "- [x] **AB**\n- [x] **CD**\n", table("built"), 0],
  ];
  const wrong = cases
    .map(([name, p, c, f, want]) => [name, judge(p, c, f, "plan.md").length, want] as const)
    .filter(([, got, want]) => got !== want)
    .map(([name, got, want]) => `"${name}": ${got} problem(s), expected ${want}`);
  const brief =
    "# Lane AB: x\n\n**Items**, a commit at each:\n\n1. **first**: a.\n2. **second**: b.\n\n**How**:\n1. **not an item**\n";
  const row = (cell: string) =>
    `| Lane or gate | Merged | Model | Items | Bug |\n|---|---|---|---|---|\n| AB | d | m | ${cell} | 0 |\n`;
  const briefCases: [string, string, string, string][] = [
    ["a row ticking every item with a proof", brief, row("1 ✓ test a; 2 ✓ output b"), ""],
    ["an item marked not landed", brief, row("1 ✓ test a; 2 ✗ not landed: followups"), ""],
    ["a proof with a semicolon", brief, row("1 ✓ a; then b; 2 ✓ c"), ""],
    ["an item missing from the row", brief, row("1 ✓ test a"), "item 2"],
    ["an item ticked without a proof", brief, row("1 ✓ test a; 2 ✓"), "item 2"],
    ["an item crossed without not landed", brief, row("1 ✓ test a; 2 ✗ later"), "item 2"],
    ["the row missing", brief, row("1 ✓ a; 2 ✓ b").replace("| AB |", "| CD |"), "no row"],
    ["a brief with no items", "# Lane AB: x\n", row("1 ✓ a"), "no numbered items"],
  ];
  for (const [name, b, f, want] of briefCases) {
    const got = judgeBrief(b, f);
    if (want ? !got.some((p) => p.includes(want)) : got.length)
      wrong.push(`"${name}": ${got.length ? got.join("; ") : "no problem"}, expected ${want || "none"}`);
    else if (want) console.log(`features-planned --self-check: planted failure named: ${got[0]}`);
  }
  if (!wrong.length)
    console.log(
      `features-planned --self-check: ${cases.length} planted plans and ${briefCases.length} planted briefs judged as expected`,
    );
  finish("features-planned --self-check", wrong);
}

if (values.has("--brief") || values.has("--findings")) {
  const briefFile =
    values.get("--brief") ?? usage("features-planned.ts --instance <dir> --brief <brief.md> --findings <findings.md>");
  const findingsFile = values.get("--findings") ?? usage("--brief needs --findings <findings.md>");
  for (const f of [briefFile, findingsFile]) if (!existsSync(join(ROOT, f))) usage(`no ${f} under ${ROOT}`);
  finish(
    "features-planned --brief",
    judgeBrief(readFileSync(join(ROOT, briefFile), "utf8"), readFileSync(join(ROOT, findingsFile), "utf8")),
  );
}

const plan = values.get("--plan") ?? usage("features-planned.ts --instance <dir> --plan <the work's plan.md>");
const FEATURES = values.get("--features") ?? `${INSTANCE}/design/features.md`;

function read(file: string): string {
  if (!existsSync(join(ROOT, file))) usage(`no ${file} under ${ROOT}`);
  return readFileSync(join(ROOT, file), "utf8");
}

const checklist = join(dirname(plan), "checklist.md");
const checklistText = existsSync(join(ROOT, checklist)) ? read(checklist) : undefined;
finish("features-planned", judge(read(plan), checklistText, read(FEATURES), plan));
