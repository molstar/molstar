// The defaults-designed check (§defaults-designed; the fourth slice's fold-back, B6): a design under `.scaffolds/seq/design/`
// whose `## Consequences` section names a layer, a layout or a part has a `## Defaults` heading (what a reader sees
// with no setting). Every design is checked except the ones in `EXEMPT`, those that failed when the check was written
// (2026-11-21): the list only shrinks, a design leaves it when it gets its Defaults section. `--worklist` prints every
// design that fails, exempt or not, so the close can see what is owed.
// Usage: defaults-designed.ts [root] [--worklist]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, trackedFiles } from "../../principles/_checks/lib.ts";

const DESIGNS = ".scaffolds/seq/design/";
/** Designs without a Defaults heading when the check was written (file names without `.md`). Only ever shrinks. */
const EXEMPT = new Set([
  "axis",
  "composition",
  "edge-cases",
  "extensibility",
  "files",
  "grids-and-rows",
  "layout-and-draw-lists",
  "layouts",
  "molecules",
  "navigators",
  "seq-data-ports",
  "testing",
  "themes",
  "track-header-parts",
]);

const argv = process.argv.slice(2);
const WORKLIST = argv.includes("--worklist");
if (WORKLIST) process.argv.splice(process.argv.indexOf("--worklist"), 1);
const { root: ROOT } = args();

/** The text of the `## <title>` section, or undefined when there is none. */
function section(text: string, title: RegExp): string | undefined {
  const lines = text.split("\n");
  const from = lines.findIndex((l) => /^## /.test(l) && title.test(l.slice(3)));
  if (from < 0) return undefined;
  let to = lines.findIndex((l, i) => i > from && /^## /.test(l));
  if (to < 0) to = lines.length;
  return lines.slice(from + 1, to).join("\n");
}

const NAMES_A_DRAWN_THING = /\b(layers?|layouts?|parts?)\b/i;
const problems: string[] = [];
const worklist: string[] = [];
for (const f of trackedFiles(ROOT)) {
  if (!f.startsWith(DESIGNS) || !f.endsWith(".md") || !existsSync(join(ROOT, f))) continue;
  const id = f.slice(DESIGNS.length, -3);
  const text = readFileSync(join(ROOT, f), "utf8");
  const consequences = section(text, /^Consequences\b/);
  if (consequences === undefined || !NAMES_A_DRAWN_THING.test(consequences)) continue;
  if (section(text, /^Defaults\b/) !== undefined) continue;
  worklist.push(f);
  if (!EXEMPT.has(id))
    problems.push(`${f}: its Consequences name a layer, a layout or a part but it has no "## Defaults" heading`);
}
if (WORKLIST)
  console.log(`defaults-designed: designs with no Defaults heading\n${worklist.map((w) => `  - ${w}`).join("\n")}`);
for (const id of EXEMPT)
  if (!worklist.includes(`${DESIGNS}${id}.md`)) problems.push(`${DESIGNS}${id}.md: passes now: drop it from EXEMPT`);
finish("defaults-designed", problems);
