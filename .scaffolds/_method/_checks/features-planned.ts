// The features-planned check (the rule "the features table moves with the merges", from Seq*'s fold-backs): a row of
// the scaffold's features table whose status says "planned" or "to come" fails once every lane of its design is merged
// in the plan (the lane that builds a feature moves its row in its own change). The rows name designs, not lanes: the
// plan's "The designs" table (`| [<design>](<file>.md) ... | <lane>, <lane>; ... |`, the method's plan template) maps
// each design to its lanes, and a lane is merged by a line `- [x] **<LANE>**`, `- [x] **<LANE>, <LANE>**` or
// `- [x] <LANE>`. A row names its designs by the links of its Design column; a row whose designs are not in the table
// cannot be checked. The status column is found by the header's index (`Feature`, `Status`, `Design`).
// The method's, not an instance's: it reads only the method's plan shape and a features table; a scaffold with no
// features table passes `--features` to another or does not run it. No instance file: the plan is always given, since
// a check never names a work folder (§execution-unreferenced).
// Usage: features-planned.ts [root] --instance .scaffolds/<name> --plan <plan.md, from the root>
//        [--features <table, from the root; default <instance>/design/features.md>]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { finish, methodArgs, usage } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, values } = methodArgs({ values: ["--plan", "--features"] });
const plan = values.get("--plan") ?? usage("features-planned.ts --instance <dir> --plan <the work's plan.md>");
const FEATURES = values.get("--features") ?? `${INSTANCE}/design/features.md`;

function read(file: string): string {
  if (!existsSync(join(ROOT, file))) usage(`no ${file} under ${ROOT}`);
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
  if (m) for (const id of (m[1] as string).match(LaneId) ?? []) merged.add(id);
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
