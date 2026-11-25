// The fold-back check (§reconcile-before-merge; README.md here). Every session record's Learned section answers the
// fold-back question on a line of its own (`Fold-back: none` or `Fold-back: proposed (foldback.md)`); a record that
// says "proposed" has a `foldback.md` in its work folder (the folder shares the record's work id); and a `foldback.md`
// in a closed work folder (its `state.md` has a `Closed` line) ends with the author's decision:
// `Decided by the author, <date>: <merged | merged in part | not merged>: "<words>"`. Records are `.<instance>/records/`
// with an instance, else `.records/`. Records written before the rule are exempted with `--allow`, since records and
// closed folders are never edited. It also notes (never fails) an open folder's `state.md` past STATE_WORDS: the state
// file holds the current state only (§state-file), and history belongs in git, the checklist and the record.
// Usage: foldback-check.ts [root] [--instance .<name>] [--allow <regex>]...
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, trackedFiles, workFolder } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, allow } = args({ allow: true });
const RECORDS = INSTANCE ? `${INSTANCE}/records/` : ".records/";
const WORK = workFolder(INSTANCE);
const ANSWER = /^Fold-back: (?:none|proposed \(foldback\.md\))\s*$/m;
const DECISION = /^\**Decided by the author\**, \d{4}-\d{2}-\d{2}: (?:merged|merged in part|not merged):/m;
const CLOSED = /^(?:#+\s*)?Closed\b/m;
const STATE_WORDS = 600;

const read = (f: string): string => {
  try {
    return readFileSync(join(ROOT, f), "utf8");
  } catch {
    return ""; // deleted in the working tree
  }
};
/** The text of a `## <name>` section, up to the next `## ` heading. */
const section = (text: string, name: string): string | undefined => {
  const m = new RegExp(`^## ${name}\\b[^\\n]*\\n([\\s\\S]*?)(?=^## |(?![\\s\\S]))`, "m").exec(text);
  return m?.[1];
};

const problems: string[] = [];
const notes: string[] = [];
const tracked = trackedFiles(ROOT);
for (const f of tracked) {
  if (allow.some((re) => re.test(f))) continue;
  if (f.startsWith(RECORDS) && f.endsWith(".md") && f.split("/").length > RECORDS.split("/").length) {
    const learned = section(read(f), "Learned");
    if (learned === undefined) {
      problems.push(`${f}: no Learned section (§permanent-record)`);
      continue;
    }
    const answer = ANSWER.exec(learned)?.[0];
    if (!answer) {
      problems.push(
        `${f}: Learned does not answer the fold-back question ("Fold-back: none" or "Fold-back: proposed (foldback.md)")`,
      );
      continue;
    }
    const workId = f.slice(RECORDS.length).split("/")[0] as string;
    if (answer.includes("proposed") && !existsSync(join(ROOT, WORK, workId, "foldback.md")))
      problems.push(`${f}: says a fold-back is proposed, but ${WORK}${workId}/foldback.md does not exist`);
  }
  if (f.startsWith(WORK) && f.endsWith("/state.md") && f.split("/").length === WORK.split("/").length + 1) {
    const text = read(f);
    const words = text.split(/\s+/).filter(Boolean).length;
    if (!CLOSED.test(text) && words > STATE_WORDS)
      notes.push(`${f}: ${words} words, over ${STATE_WORDS}: rewrite it to the current state only (§state-file)`);
  }
  if (f.startsWith(WORK) && f.endsWith("/foldback.md")) {
    const dir = f.slice(0, -"foldback.md".length);
    if (CLOSED.test(read(`${dir}state.md`)) && !DECISION.test(read(f)))
      problems.push(
        `${f}: the folder is closed, but the author's decision line is missing ("Decided by the author, <date>: <merged | merged in part | not merged>: ...")`,
      );
  }
}
if (notes.length) console.log(`foldback-check: ${notes.length} note(s)\n${notes.map((n) => `  - ${n}`).join("\n")}`);
finish("foldback-check", problems);
