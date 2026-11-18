// The closed-work check (§reconcile-before-merge; the manifest, "The fold-back is a step of every close"): a work
// folder under `.seq/execution/` whose `state.md` ends with a `Closed` line has a `foldback.md`, and that file records
// the author's approval on a line of its own starting `Approved by the author`. Folders closed before the rule (listed
// below) are exempt: closed folders are never edited.
// Usage: closed-work.ts [root]
import { existsSync, readdirSync, readFileSync, statSync } from "node:fs";
import { join } from "node:path";
import { args, finish } from "../oz/_checks/lib.ts";

const WORK = ".seq/execution";
/** Closed before the rule was written (2026-11-18): their fold-backs, if any, were approved in the session record. */
const BEFORE_THE_RULE = new Set(["adopt-principles", "first-slice"]);

const { root: ROOT } = args();
const problems: string[] = [];
const work = join(ROOT, WORK);
for (const id of existsSync(work) ? readdirSync(work) : []) {
  const dir = join(work, id);
  if (!statSync(dir).isDirectory() || BEFORE_THE_RULE.has(id)) continue;
  const state = join(dir, "state.md");
  if (!existsSync(state)) continue;
  const lines = readFileSync(state, "utf8").split("\n").map((l) => l.trim()).filter(Boolean);
  if (!/^Closed\b/.test(lines[lines.length - 1] ?? "")) continue; // still open
  const foldback = join(dir, "foldback.md");
  if (!existsSync(foldback)) {
    problems.push(`${WORK}/${id}: closed without a foldback.md`);
    continue;
  }
  if (!/^\**Approved by the author\b/m.test(readFileSync(foldback, "utf8")))
    problems.push(`${WORK}/${id}/foldback.md: no "Approved by the author" line`);
}
finish("closed-work", problems);
