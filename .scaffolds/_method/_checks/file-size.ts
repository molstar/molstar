// The module-size check (§module-size): a tracked source file in the scaffold's checked folders has at most `hard`
// lines (tests `testHard`), unless the baseline lists it, and then at most its listed count: the baseline only
// shrinks. A file past `soft` is printed as a note (a question at review), not a failure. `--update` lowers the
// baseline to what the files now hold and drops the files back under the bound; it never raises a count or adds a
// file (the orchestrator runs it after a split; a lane never edits the baseline, §lanes-dont-edit-checks).
// The method's, not an instance's: the bound is the package's rule; which folders it governs is the scaffold's. The
// instance's part: `<instance>/_checks/file-size.json`, required, `{ "checked": "<regex over tracked paths>" }` and
// optionally `soft`, `hard`, `testHard` (default 800, 1500, 2500, §module-size's); and
// `<instance>/_checks/file-size-baseline.json`, the files past the bound when the check was written (absent, none).
// Usage: file-size.ts [root] --instance .scaffolds/<name> [--update]
import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { finish, instanceFile, methodArgs, trackedFiles, usage } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, flags } = methodArgs({ flags: ["--update"] });
const update = flags.has("--update");
const config = instanceFile<{ checked?: string; soft?: number; hard?: number; testHard?: number }>(
  ROOT,
  INSTANCE,
  "file-size.json",
);
if (!config.checked) usage(`${INSTANCE}/_checks/file-size.json: "checked" (a regex over tracked paths) is required`);
/** The folders checked: the scaffold's packages. */
const CHECKED = new RegExp(config.checked);
const SKIPPED = /(^|\/)(node_modules|dist|lib|build)\//;
const TEST = /(\/_test\/|\.test\.|\.quick\.|\.heap\.)/;
const SOFT = config.soft ?? 800;
const HARD = config.hard ?? 1500;
const TEST_HARD = config.testHard ?? 2500;
const BASELINE = `${INSTANCE}/_checks/file-size-baseline.json`;

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
