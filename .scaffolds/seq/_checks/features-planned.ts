// The features-planned check (the fourth slice's fold-back, B7): a row of `.seq/design/features.md` whose status says
// "planned" or "to come" fails once every lane of its design is merged in the plan (§features-move-with-merges: the
// lane that builds a feature moves its row in its own change). The rows name designs, not lanes: the plan's "The
// designs" table (`| [<design>](<file>.md) ... | <lane>, <lane>; ... |`) maps each design to its lanes, and a lane is
// merged by a line `- [x] **<LANE>**`, `- [x] **<LANE>, <LANE>**` or `- [x] <LANE>`. A row names its designs by the
// links of its Design column; a row whose designs are not in the table cannot be checked. The status column is found
// by the header's index. The plan is the fifth slice's by default.
// Usage: features-planned.ts [root] [--plan <path, relative to the root>]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish } from "../oz/_checks/lib.ts";

const FEATURES = ".seq/design/features.md";
const argv = process.argv.slice(2);
let plan: string | undefined; // the slice's plan, always given: a check never names a work folder (§execution-unreferenced)
const at = argv.indexOf("--plan");
if (at >= 0) {
  const v = argv[at + 1];
  if (v === undefined || v.startsWith("--")) {
    console.error("usage: --plan needs a value");
    process.exit(2);
  }
  plan = v;
  argv.splice(at, 2);
}
if (plan === undefined) {
  console.error("usage: features-planned.ts --plan <the slice's plan.md>");
  process.exit(2);
}
process.argv.splice(2, process.argv.length, ...argv);
const { root: ROOT } = args();

function read(file: string): string {
  if (!existsSync(join(ROOT, file))) {
    console.error(`usage: no ${file} under ${ROOT}`);
    process.exit(2);
  }
  return readFileSync(join(ROOT, file), "utf8");
}

const LaneId = /\b[A-Z][A-Z0-9]+\b/g;
const planText = read(plan);

/** The lanes each design of the plan's table has. */
const lanesOf = new Map<string, string[]>();
for (const line of planText.split("\n")) {
  if (!line.startsWith("|")) continue;
  const cells = line.split("|").slice(1, -1).map((c) => c.trim());
  const design = /^\[([^\]]+)\]\(/.exec(cells[0] ?? "")?.[1];
  if (design && cells[1]) lanesOf.set(design, cells[1].match(LaneId) ?? []);
}

/** The lanes the plan marks merged. */
const merged = new Set<string>();
for (const line of planText.split("\n")) {
  const m = /^\s*- \[x\]\s+\**((?:[A-Z][A-Z0-9]*)(?:\s*[,;]\s*[A-Z][A-Z0-9]*)*)/.exec(line);
  if (m) for (const id of m[1].match(LaneId) ?? []) merged.add(id);
}

const problems: string[] = [];
let header: string[] = [];
for (const line of read(FEATURES).split("\n")) {
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
    problems.push(`${cells[0]}: status "${status.slice(0, 40)}" but every lane of ${designs.join(", ")} (${lanes.join(", ")}) is merged in ${plan}`);
}
finish("features-planned", problems);
