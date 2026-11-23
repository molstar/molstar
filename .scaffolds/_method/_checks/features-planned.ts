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
// Usage: features-planned.ts [root] --instance .scaffolds/<name> --plan <plan.md, from the root>
//        [--features <table, from the root; default <instance>/design/features.md>] | [--self-check]
//        (self-check: planted plans, checklists and tables, each failing or passing as it should)
import { existsSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { finish, methodArgs, usage } from "./lib.ts";

const {
  root: ROOT,
  instance: INSTANCE,
  values,
  flags,
} = methodArgs({ values: ["--plan", "--features"], flags: ["--self-check"] });

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
  if (!wrong.length) console.log(`features-planned --self-check: ${cases.length} planted plans judged as expected`);
  finish("features-planned --self-check", wrong);
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
