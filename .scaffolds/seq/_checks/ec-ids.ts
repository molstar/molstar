// The edge-case id check (§mechanical-drift-checks; the manifest, "Folded back from the first slice"): every `EC-<n>`
// cited in a tracked file under `packages/seq` (built output in `lib/` aside) names an item of the catalog,
// `.seq/sources/seqstar/edge-cases.md`, where an item is a top-level numbered line (`143. ...`). A range written
// `EC-<a>–<b>` or `EC-<a> to EC-<b>` cites each id in it. The ids are numbers only, so code citing them names no
// `.seq/` path; the design's `edge-cases.md` only refers to the catalog.
// Coverage (the third slice's fold-back, B4): a test names the catalog item it tests (`(EC-43)` in its title or a
// comment of its file under `_test/`). A feature marked `done` in `.seq/design/features.md` covers its Cases: each
// must be named by a test or have an outcome in `.seq/design/edge-cases.md` ("Outcomes of ..."), else it fails,
// unless the baseline lists it (rows marked done before the rule; the baseline only shrinks, and a lane never edits
// it, §lanes-dont-edit-checks). The items a built row (`done` or `building`) covers that no test names are printed as
// the close's worklist (`--worklist`).
// Usage: ec-ids.ts [root] [--worklist]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, trackedFiles } from "../oz/_checks/lib.ts";

const CATALOG = ".seq/sources/seqstar/edge-cases.md";
const SCOPE = "packages/seq/";

const argv = process.argv.slice(2);
const WORKLIST = argv.includes("--worklist");
if (WORKLIST) process.argv.splice(process.argv.indexOf("--worklist"), 1);
const { root: ROOT } = args();
const FEATURES = ".seq/design/features.md";
const OUTCOMES = ".seq/design/edge-cases.md";
const BASELINE = ".seq/_checks/ec-coverage-baseline.json";
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

/** `EC-12`, `EC-9–14` / `EC-9-14`, `EC-9 to EC-14`: the ids each citation names. */
const CITE = /\bEC-(\d+)(?:\s*[–-]\s*(\d+)|\s+to\s+EC-(\d+))?/g;

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
