// The closed-work check (§reconcile-before-merge, §plan-and-checklist, §state-file; the rule "the fold-back is a step of
// every close"): a work folder under `<instance>/execution/` whose `state.md` ends with a `Closed` line
// - has a `foldback.md` that records the author's decision on a line of its own:
//   `Decided by the author, <date>: <merged | merged in part | not merged>:`. Writing it is required, merging it is not
//   (an exploration's lessons stay in its archive);
// - has no unticked line (`- [ ]`) in its `checklist.md`: Finalize runs, cuts or carries each item and ticks it so.
//   A folder closed before checklists existed has none and passes as before.
// An open folder's `state.md` past STATE_WORDS words is noted, never failed: it holds the current state only, and
// history belongs in git, the checklist's ticks, `findings.md` and the record.
// The method's, not an instance's: every scaffold closes its work the same way. The instance's part is
// `<instance>/_checks/closed-work.json`, `{ "exempt": ["<work id>", ...] }`: folders closed before the rule, never
// edited since (closed folders are not); absent, none is exempt.
// Usage: closed-work.ts [root] --instance .scaffolds/<name> [--self-check] (self-check: planted folders, one per rule,
// each failing or noted as it should, and the folders that must pass passing).
import { existsSync, readdirSync, readFileSync, statSync } from "node:fs";
import { join } from "node:path";
import { finish, instanceFile, methodArgs } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, flags } = methodArgs({ flags: ["--self-check"] });
const WORK = `${INSTANCE}/execution`;
const STATE_WORDS = 600;

/** The decision line: who, when, and one of the three outcomes. */
const DECISION = /^\**Decided by the author\**, \d{4}-\d{2}-\d{2}: (?:merged|merged in part|not merged):/m;
/** An unticked checklist line. */
const UNTICKED = /^\s*- \[ \]\s*(.*)$/gm;

/** A work folder's files, as text; a missing file is undefined. */
interface Folder {
  state?: string;
  foldback?: string;
  checklist?: string;
}

/** The problems and notes of one work folder. */
function judge(id: string, f: Folder): { problems: string[]; notes: string[] } {
  const problems: string[] = [];
  const notes: string[] = [];
  if (f.state === undefined) return { problems, notes };
  const lines = f.state.split("\n").map((l) => l.trim()).filter(Boolean);
  if (!/^Closed\b/.test(lines[lines.length - 1] ?? "")) {
    const words = f.state.split(/\s+/).filter(Boolean).length;
    if (words > STATE_WORDS)
      notes.push(`${WORK}/${id}/state.md: ${words} words, over ${STATE_WORDS}: rewrite it to the current state only (§state-file)`);
    return { problems, notes };
  }
  if (f.foldback === undefined) problems.push(`${WORK}/${id}: closed without a foldback.md`);
  else if (!DECISION.test(f.foldback))
    problems.push(`${WORK}/${id}/foldback.md: no "Decided by the author, <date>: <merged | merged in part | not merged>:" line`);
  for (const m of (f.checklist ?? "").matchAll(UNTICKED))
    problems.push(`${WORK}/${id}/checklist.md: closed with an unticked item "${(m[1] as string).slice(0, 60)}": run it, cut it by the author's word or carry it, and tick it so`);
  return { problems, notes };
}

if (flags.has("--self-check")) {
  const closed = "Done.\n\nClosed 2026-01-01: records/w/2026-01-01-w.md\n";
  const decided = 'Decided by the author, 2026-01-01: merged: "yes".\n';
  const long = `# State\n\n${"word ".repeat(STATE_WORDS + 1)}\n`;
  const cases: [string, Folder, number, number][] = [
    ["closed, decided, all ticked", { state: closed, foldback: decided, checklist: "- [x] A (abc)\n- [x] FINALIZE (def)\n" }, 0, 0],
    ["closed before checklists", { state: closed, foldback: decided }, 0, 0],
    ["closed, an unticked item", { state: closed, foldback: decided, checklist: "- [x] A (abc)\n- [ ] GATE 2\n" }, 1, 0],
    ["closed, FINALIZE unticked", { state: closed, foldback: decided, checklist: "- [ ] FINALIZE: the close\n" }, 1, 0],
    ["closed, no foldback", { state: closed, checklist: "- [x] A\n" }, 1, 0],
    ["closed, no decision", { state: closed, foldback: "# Fold-back\n" }, 1, 0],
    ["open, unticked, short", { state: "# State\n\nNext: A.\n", checklist: "- [ ] A\n" }, 0, 0],
    ["open, long state", { state: long, checklist: "- [ ] A\n" }, 0, 1],
    ["closed, long state", { state: `${long}${closed}`, foldback: decided }, 0, 0],
  ];
  const wrong: string[] = [];
  for (const [name, folder, p, n] of cases) {
    const r = judge("planted", folder);
    if (r.problems.length !== p || r.notes.length !== n)
      wrong.push(`"${name}": ${r.problems.length} problem(s) and ${r.notes.length} note(s), expected ${p} and ${n}`);
  }
  if (!wrong.length) console.log(`closed-work --self-check: ${cases.length} planted folders judged as expected`);
  finish("closed-work --self-check", wrong);
}

const { exempt = [] } = instanceFile<{ exempt?: string[] }>(ROOT, INSTANCE, "closed-work.json", {});
const BEFORE_THE_RULE = new Set(exempt);
const read = (dir: string, file: string): string | undefined =>
  existsSync(join(dir, file)) ? readFileSync(join(dir, file), "utf8") : undefined;

const problems: string[] = [];
const notes: string[] = [];
const work = join(ROOT, WORK);
for (const id of existsSync(work) ? readdirSync(work) : []) {
  const dir = join(work, id);
  if (!statSync(dir).isDirectory() || BEFORE_THE_RULE.has(id)) continue;
  const r = judge(id, { state: read(dir, "state.md"), foldback: read(dir, "foldback.md"), checklist: read(dir, "checklist.md") });
  problems.push(...r.problems);
  notes.push(...r.notes);
}
for (const id of BEFORE_THE_RULE)
  if (!existsSync(join(work, id))) problems.push(`closed-work.json: exempt "${id}" has no folder in ${WORK}: drop it`);
if (notes.length) console.log(`closed-work: ${notes.length} note(s)\n${notes.map((n) => `  - ${n}`).join("\n")}`);
finish("closed-work", problems);
