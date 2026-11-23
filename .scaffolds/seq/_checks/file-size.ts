// The module-size check (§module-size): a tracked source file in the checked folders has at most HARD lines (tests
// TEST_HARD), unless the baseline lists it, and then at most its listed count: the baseline only shrinks. A file past
// SOFT is printed as a note (a question at review), not a failure. `--update` lowers the baseline to what the files now
// hold and drops the files back under the bound; it never raises a count or adds a file (the orchestrator runs it after
// a split; a lane never edits the baseline, §lanes-dont-edit-checks).
// Usage: file-size.ts [root] [--update]
import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { finish, trackedFiles } from "../oz/_checks/lib.ts";

/** The folders checked: Seq*'s packages. */
const CHECKED = /^packages\/seq\/.*\.(ts|tsx|mts|mjs|js)$/;
const SKIPPED = /(^|\/)(node_modules|dist|lib|build)\//;
const TEST = /(\/_test\/|\.test\.|\.quick\.|\.heap\.)/;
const SOFT = 800;
const HARD = 1500;
const TEST_HARD = 2500;
const BASELINE = ".seq/_checks/file-size-baseline.json";

const argv = process.argv.slice(2);
const update = argv.includes("--update");
const rest = argv.filter((a) => a !== "--update");
if (rest.some((a) => a.startsWith("-")) || rest.length > 1) {
  console.error("usage: file-size.ts [root] [--update]");
  process.exit(2);
}
const ROOT = rest[0] ?? ".";

const lineCount = (text: string): number => {
  if (!text) return 0;
  const n = text.split("\n").length;
  return text.endsWith("\n") ? n - 1 : n;
};

const baselinePath = join(ROOT, BASELINE);
const baseline: Record<string, number> = existsSync(baselinePath)
  ? JSON.parse(readFileSync(baselinePath, "utf8"))
  : {};

const problems: string[] = [];
const notes: string[] = [];
const next: Record<string, number> = {};
const seen = new Set<string>();
for (const f of trackedFiles(ROOT)) {
  if (!CHECKED.test(f) || SKIPPED.test(f)) continue;
  let text: string;
  try {
    text = readFileSync(join(ROOT, f), "utf8");
  } catch {
    continue; // deleted in the working tree
  }
  seen.add(f);
  const lines = lineCount(text);
  const hard = TEST.test(f) ? TEST_HARD : HARD;
  const listed = baseline[f];
  if (listed !== undefined) {
    if (lines > listed)
      problems.push(`${f}: ${lines} lines, grown past its baseline ${listed} (split it; the baseline only shrinks)`);
    if (lines > hard) next[f] = Math.min(lines, listed);
    else if (!update) notes.push(`${f}: ${lines} lines, back under ${hard}: \`--update\` drops it from the baseline`);
    if (lines < listed && lines > hard && !update)
      notes.push(`${f}: ${lines} lines, under its baseline ${listed}: \`--update\` lowers it`);
  } else if (lines > hard) {
    problems.push(`${f}: ${lines} lines, over ${hard} (§module-size: split it by what changes together)`);
  } else if (lines > SOFT && !TEST.test(f)) {
    notes.push(`${f}: ${lines} lines, over ${SOFT}: say why it stays whole, or split it`);
  }
}
for (const f of Object.keys(baseline))
  if (!seen.has(f) && !update) notes.push(`${f}: listed in the baseline, no longer tracked: \`--update\` drops it`);

if (update) {
  const sorted = Object.fromEntries(Object.entries(next).sort(([a], [b]) => a.localeCompare(b)));
  writeFileSync(baselinePath, `${JSON.stringify(sorted, null, 2)}\n`);
  console.log(`file-size: baseline written, ${Object.keys(sorted).length} file(s)`);
}
if (notes.length) console.log(`file-size: ${notes.length} note(s)\n${notes.map((n) => `  - ${n}`).join("\n")}`);
finish("file-size", problems);
