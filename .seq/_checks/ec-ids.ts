// The edge-case id check (§mechanical-drift-checks; the manifest, "Folded back from the first slice"): every `EC-<n>`
// cited in a tracked file under `packages/seq` (built output in `lib/` aside) names an item of the catalog,
// `.seq/sources/seqstar/edge-cases.md`, where an item is a top-level numbered line (`143. ...`). A range written
// `EC-<a>–<b>` or `EC-<a> to EC-<b>` cites each id in it. The ids are numbers only, so code citing them names no
// `.seq/` path; the design's `edge-cases.md` only refers to the catalog.
// Usage: ec-ids.ts [root]
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, trackedFiles } from "../oz/_checks/lib.ts";

const CATALOG = ".seq/sources/seqstar/edge-cases.md";
const SCOPE = "packages/seq/";

const { root: ROOT } = args();

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
    }
  }
}
if (!problems.length) console.log(`ec-ids: ${cited} citation(s) against ${known.size} catalog items`);
finish("ec-ids", problems);
