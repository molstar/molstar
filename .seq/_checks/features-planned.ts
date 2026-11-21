// The features-planned check (the fourth slice's fold-back, B7): a row of `.seq/design/features.md` whose status says
// "planned" or "to come" while the row names a lane that the plan marks merged (`- [x] **<LANE>**`) fails: the lane
// that builds a feature moves its row in its own change (§features-move-with-merges). A row names a lane by its id
// as a word (`GG`, `GC1`). The plan is the fifth slice's by default.
// Usage: features-planned.ts [root] [--plan <path, relative to the root>]
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish } from "../oz/_checks/lib.ts";

const FEATURES = ".seq/design/features.md";
const argv = process.argv.slice(2);
let plan = ".seq/execution/fifth-slice/plan.md";
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
process.argv.splice(2, process.argv.length, ...argv);
const { root: ROOT } = args();

const merged = new Set<string>();
for (const m of readFileSync(join(ROOT, plan), "utf8").matchAll(/^\s*- \[x\] \*\*([A-Z][A-Z0-9]*)\*\*/gm))
  merged.add(m[1] as string);

const problems: string[] = [];
for (const line of readFileSync(join(ROOT, FEATURES), "utf8").split("\n")) {
  if (!line.startsWith("|")) continue;
  const cells = line
    .split("|")
    .slice(1, -1)
    .map((c) => c.trim());
  const status = cells.find((c) => /^(planned|building|done|deferred)\b/.test(c));
  if (status === undefined || !/\b(planned|to come)\b/.test(status)) continue;
  for (const lane of merged)
    if (new RegExp(`\\b${lane}\\b`).test(line))
      problems.push(`${cells[0]}: status "${status.slice(0, 40)}" but ${lane} is merged in ${plan}`);
}
finish("features-planned", problems);
