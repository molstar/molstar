// The closed-work check (§reconcile-before-merge; the rule "the fold-back is a step of every close"): a work folder
// under `<instance>/execution/` whose `state.md` ends with a `Closed` line has a `foldback.md`, and that file records
// the author's decision on a line of its own: `Decided by the author, <date>: <merged | merged in part | not merged>:`.
// Writing it is required, merging it is not (an exploration's lessons stay in its archive).
// The method's, not an instance's: every scaffold closes its work the same way. The instance's part is
// `<instance>/_checks/closed-work.json`, `{ "exempt": ["<work id>", ...] }`: folders closed before the rule, never
// edited since (closed folders are not); absent, none is exempt.
// Usage: closed-work.ts [root] --instance .scaffolds/<name>
import { existsSync, readdirSync, readFileSync, statSync } from "node:fs";
import { join } from "node:path";
import { finish, instanceFile, methodArgs } from "./lib.ts";

const { root: ROOT, instance: INSTANCE } = methodArgs();
const WORK = `${INSTANCE}/execution`;
const { exempt = [] } = instanceFile<{ exempt?: string[] }>(ROOT, INSTANCE, "closed-work.json", {});
const BEFORE_THE_RULE = new Set(exempt);

/** The decision line: who, when, and one of the three outcomes. */
const DECISION = /^\**Decided by the author\**, \d{4}-\d{2}-\d{2}: (?:merged|merged in part|not merged):/m;

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
  if (!DECISION.test(readFileSync(foldback, "utf8")))
    problems.push(`${WORK}/${id}/foldback.md: no "Decided by the author, <date>: <merged | merged in part | not merged>:" line`);
}
for (const id of BEFORE_THE_RULE)
  if (!existsSync(join(work, id))) problems.push(`closed-work.json: exempt "${id}" has no folder in ${WORK}: drop it`);
finish("closed-work", problems);
