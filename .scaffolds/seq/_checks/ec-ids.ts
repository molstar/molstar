// The edge-case id check (§mechanical-drift-checks; the manifest's §ec-tested-by-name): every `EC-<n>`
// cited in a tracked file under `packages/seq` (built output in `lib/` aside) names an item of the catalog,
// `.scaffolds/seq/sources/seqstar/edge-cases.md`, where an item is a top-level numbered line (`143. ...`). A range written
// `EC-<a>–<b>` or `EC-<a> to EC-<b>` cites each id in it. The ids are numbers only, so code citing them names no
// `.scaffolds/seq/` path; the design's `edge-cases.md` only refers to the catalog.
// Coverage (the third slice's fold-back, B4): a test names the catalog item it tests (`(EC-43)` in its title or a
// comment of its file under `_test/`). A feature marked `done` in `.scaffolds/seq/design/features.md` covers its Cases: each
// must be named by a test or have an outcome in `.scaffolds/seq/design/edge-cases.md` ("Outcomes of ..."), else it fails,
// unless the baseline lists it (rows marked done before the rule; the baseline only shrinks, and a lane never edits
// it, §lanes-dont-edit-checks). The items a built row (`done` or `building`) covers that no test names are printed as
// the close's worklist (`--worklist`).
// A rewrite keeps its ids (the sixth slice's fold-back, B2): `--compare <ref>` (a mode of its own) fails an `EC-<n>`
// cited in a file of the scope at the git ref `<ref>` (a lane's base) and no longer cited in that file, or in the file
// it was renamed to (`git diff -M` against the working tree); a deleted file's ids are not asked for. Per file, because
// a rewrite dropped EC-169 from one comment while 48 other files still cited it. `--self-check` runs planted cases.
// Usage: ec-ids.ts [root] [--worklist] | ec-ids.ts [root] --compare <ref> | ec-ids.ts --self-check
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, git, trackedFiles } from "../../principles/_checks/lib.ts";

const CATALOG = ".scaffolds/seq/sources/seqstar/edge-cases.md";
const SCOPE = "packages/seq/";

/** Takes a boolean flag, or a flag and its value, off the arguments before the package's parser sees them. */
function take(flag: string, valued = false): string | boolean {
  const i = process.argv.indexOf(flag, 2);
  if (i < 0) return false;
  const v = valued ? process.argv[i + 1] : true;
  if (valued && (v === undefined || v === "" || (v as string).startsWith("--"))) {
    console.error(`usage: ${flag} needs a value`);
    process.exit(2);
  }
  process.argv.splice(i, valued ? 2 : 1);
  return v as string | boolean;
}
const WORKLIST = take("--worklist") === true;
const SELF = take("--self-check") === true;
const COMPARE = take("--compare", true);
const { root: ROOT } = args();

/** `EC-12`, `EC-9–14` / `EC-9-14`, `EC-9 to EC-14`: the ids each citation names. */
const CITE = /\bEC-(\d+)(?:\s*[–-]\s*(\d+)|\s+to\s+EC-(\d+))?/g;

/** The ids a text cites, a range expanded (a backwards range names none). */
function citedIn(text: string): Set<number> {
  const out = new Set<number>();
  for (const m of text.matchAll(CITE))
    for (let id = Number(m[1]); id <= Number(m[2] ?? m[3] ?? m[1]); id++) out.add(id);
  return out;
}
/** A file of the scope changed since the ref: its path then and now (`now` absent when it was deleted). */
interface Changed {
  then: string;
  now?: string;
  before: string;
  after: string;
}
/** Each id cited in a changed file at the ref and no longer cited in it (or its renamed successor), as a problem. */
function droppedPerFile(changes: readonly Changed[], ref: string): string[] {
  const out: string[] = [];
  for (const c of changes) {
    if (c.now === undefined) continue; // deleted: its ids go with it
    const after = citedIn(c.after);
    for (const id of [...citedIn(c.before)].sort((a, b) => a - b))
      if (!after.has(id))
        out.push(`${c.now}: cited EC-${id} at ${ref}${c.now === c.then ? "" : ` (as ${c.then})`} and no longer: a rewrite dropped it, or say why it went`);
  }
  return out;
}

if (SELF) {
  const wrong: string[] = [];
  const same = (got: unknown[], want: unknown[], what: string) => {
    if (got.join("\n") !== want.join("\n")) wrong.push(`${what}: got [${got}], expected [${want}]`);
  };
  same([...citedIn("(EC-3) and EC-9–11, EC-20 to EC-21; EC-5-4")], [3, 9, 10, 11, 20, 21], "citations read");
  // the sixth slice's S1E: EC-169 dropped from one comment of rows-levels.ts while other files still cite it
  const rows = "packages/seq/core/src/render/rows-levels.ts";
  const s1e: Changed[] = [
    { then: rows, now: rows, before: "// tiles at every level (EC-169, EC-12)", after: "// tiles at every level (Render, Rows) EC-12" },
    { then: "packages/seq/core/src/a.ts", now: "packages/seq/core/src/a.ts", before: "EC-169", after: "EC-169 still" },
  ];
  same(droppedPerFile(s1e, "base").map((p) => p.split(":")[0]), [rows], "an id dropped from one file, still cited in another");
  same(droppedPerFile([{ then: "x.ts", now: "y.ts", before: "EC-7", after: "EC-7" }], "base"), [], "a rename keeping its ids");
  same(droppedPerFile([{ then: "x.ts", now: "y.ts", before: "EC-7", after: "" }], "base").length ? [1] : [], [1], "a rename dropping one");
  same(droppedPerFile([{ then: "x.ts", before: "EC-7", after: "" }], "base"), [], "a deleted file");
  same(droppedPerFile([{ then: "x.ts", now: "x.ts", before: "EC-1–3", after: "EC-1, EC-3" }], "base").length ? [1] : [], [1], "an id dropped from a range");
  if (!wrong.length) console.log("ec-ids --self-check: 6 planted cases as expected");
  finish("ec-ids --self-check", wrong);
}

if (typeof COMPARE === "string") {
  const tryGit = (...a: string[]): string => {
    try {
      return git(ROOT, ...a);
    } catch (e) {
      console.error(`usage: --compare ${COMPARE}: ${(e as Error).message.split("\n")[0]}`);
      process.exit(2);
    }
  };
  const fields = tryGit("diff", "-M", "--name-status", "-z", COMPARE, "--", SCOPE).split("\0");
  const changes: Changed[] = [];
  const inScope = (f: string) => f.startsWith(SCOPE) && !/(^|\/)lib\//.test(f);
  for (let i = 0; i < fields.length - 1; ) {
    const status = fields[i++] as string;
    const then = fields[i++] as string;
    const now = /^[RC]/.test(status) ? (fields[i++] as string) : status === "D" ? undefined : then;
    if (!/^[MRD]/.test(status) || !inScope(then)) continue;
    const before = tryGit("show", `${COMPARE}:${then}`);
    const after = now === undefined || !existsSync(join(ROOT, now)) ? "" : readFileSync(join(ROOT, now), "utf8");
    if (before.includes("\0") || after.includes("\0")) continue; // binary
    changes.push({ then, now, before, after });
  }
  const problems = droppedPerFile(changes, COMPARE);
  if (!problems.length) console.log(`ec-ids --compare: ${changes.length} changed file(s) under ${SCOPE} since ${COMPARE} keep their ids`);
  finish("ec-ids --compare", problems);
}
const FEATURES = ".scaffolds/seq/design/features.md";
const OUTCOMES = ".scaffolds/seq/design/edge-cases.md";
const BASELINE = ".scaffolds/seq/_checks/ec-coverage-baseline.json";
const TEST_FILE = /(\/_test\/|\.test\.|\.quick\.)/;

let catalog: string;
try {
  catalog = readFileSync(join(ROOT, CATALOG), "utf8");
} catch {
  console.error(`usage: no catalog at ${CATALOG} under ${ROOT}`);
  process.exit(2);
}
const known = new Set<number>();
for (const m of catalog.matchAll(/^(\d+)\. /gm)) known.add(Number(m[1]));
if (known.size === 0) {
  console.error(`usage: ${CATALOG} lists no numbered items`);
  process.exit(2);
}

const problems: string[] = [];
let cited = 0;
const tested = new Set<number>();
for (const f of trackedFiles(ROOT)) {
  if (!f.startsWith(SCOPE) || /(^|\/)lib\//.test(f)) continue;
  let text: string;
  try {
    text = readFileSync(join(ROOT, f), "utf8");
  } catch {
    continue; // deleted in the working tree
  }
  if (text.includes("\0")) continue; // binary
  for (const m of text.matchAll(CITE)) {
    const from = Number(m[1]);
    const to = Number(m[2] ?? m[3] ?? m[1]);
    const line = text.slice(0, m.index).split("\n").length;
    if (to < from) {
      problems.push(`${f}:${line} cites a backwards range ${m[0]}`);
      continue;
    }
    for (let id = from; id <= to; id++) {
      cited++;
      if (!known.has(id)) problems.push(`${f}:${line} cites EC-${id}, which ${CATALOG} does not list`);
      else if (TEST_FILE.test(f)) tested.add(id);
    }
  }
}

/** `1–8`, `68–71, 77`, `—`: the ids a Cases cell or an outcome row names. */
const idsOf = (cell: string): number[] => {
  const out: number[] = [];
  for (const m of cell.replace(/\([^)]*\)/g, "").matchAll(/(\d+)(?:\s*[–-]\s*(\d+))?/g))
    for (let id = Number(m[1]); id <= Number(m[2] ?? m[1]); id++) out.push(id);
  return out;
};
const read = (rel: string): string => readFileSync(join(ROOT, rel), "utf8");

/** Items with an outcome: the first column of each "Kept" table and each "Changed" bullet, in every slice's section. */
const outcomes = new Set<number>();
for (const section of read(OUTCOMES).split(/^## Outcomes/m).slice(1)) {
  const [kept, rest = ""] = section.split("**Changed**");
  for (const line of kept.split("\n"))
    if (/^\| [\d ,–-]/.test(line)) for (const id of idsOf(line.split("|")[1])) outcomes.add(id);
  for (const m of rest.split("**Not met")[0].matchAll(/^- EC-(\d+)(?:, EC-(\d+))?/gm)) {
    outcomes.add(Number(m[1]));
    if (m[2]) outcomes.add(Number(m[2]));
  }
}

const baseline: Record<string, number[]> = existsSync(join(ROOT, BASELINE)) ? JSON.parse(read(BASELINE)) : {};
const notes: string[] = [];
const worklist: string[] = [];
let header: string[] = [];
for (const line of read(FEATURES).split("\n")) {
  if (!line.startsWith("| ") || line.startsWith("|---")) continue;
  const cells = line.split("|").slice(1, -1).map((c) => c.trim());
  if (cells[0] === "Feature") {
    header = cells;
    continue;
  }
  const status = cells[header.indexOf("Status")] ?? "";
  const cases = header.includes("Cases") ? idsOf(cells[header.indexOf("Cases")] ?? "") : [];
  if (!/^(done|building)/.test(status)) continue;
  const untested = cases.filter((id) => !tested.has(id));
  if (untested.length) worklist.push(`${cells[0]}: EC-${untested.join(", EC-")}`);
  if (!status.startsWith("done")) continue;
  const open = cases.filter((id) => !tested.has(id) && !outcomes.has(id));
  const listed = new Set(baseline[cells[0]] ?? []);
  for (const id of open)
    if (!listed.has(id))
      problems.push(`${FEATURES}: "${cells[0]}" is done, but no test names EC-${id} and it has no outcome in ${OUTCOMES}`);
  for (const id of listed)
    if (!open.includes(id)) notes.push(`${BASELINE}: "${cells[0]}" EC-${id} is covered now: drop it from the baseline`);
}

if (!problems.length) console.log(`ec-ids: ${cited} citation(s) against ${known.size} catalog items; ${tested.size} named by a test`);
if (notes.length) console.log(`ec-ids: ${notes.length} note(s)\n${notes.map((n) => `  - ${n}`).join("\n")}`);
if (WORKLIST) console.log(`ec-ids: built rows' items no test names\n${worklist.map((w) => `  - ${w}`).join("\n")}`);
finish("ec-ids", problems);
