// The no-design-citations check (stage 1b of the sixth slice): this folder is temporary, so nothing under `packages/`
// (or the command-line packages in `cli/`) names a design or a source file of the scaffold (`.seq/design/*.md`,
// `.seq/sources/**/*.md`, read at run time) or a path into `.seq/`; a decision is cited by the spec section that holds
// it. A file name that also names a tracked file outside the scaffold (a package's own `_spec/README.md`, `format.md`
// of another package's spec) is that package's own and passes. The one pointer to the scaffold is
// `packages/seq/AGENTS.md`, until the scaffold is deleted.
// Usage: no-design-citations.ts [root]
import { existsSync, readFileSync } from "node:fs";
import { basename, join } from "node:path";
import { args, finish, trackedFiles } from "../oz/_checks/lib.ts";

const { root: ROOT } = args();
const POINTER = "packages/seq/AGENTS.md";
/** Files that name the scaffold when the check was written (2026-11-23) and that a lane may not edit (the
 * orchestrator's hot files): the list only shrinks, a file leaves it when it is clean (the check then says so). */
const PENDING = new Set(["packages/seq/TODO.md", "packages/seq/_measure/sizes.mjs"]);
const TEXT = /\.(ts|tsx|mts|mjs|js|json|md|mdx|css|scss|html|yml|yaml|txt)$/;

const files = trackedFiles(ROOT).filter((f) => existsSync(join(ROOT, f)));

/** The scaffold's design and source file names (without the folder, README.md left out: every folder has one). */
const scaffold = new Set<string>();
for (const f of files) {
  if (!f.endsWith(".md") || basename(f) === "README.md") continue;
  if (f.startsWith(".seq/design/") || f.startsWith(".seq/sources/")) scaffold.add(basename(f));
}
/** Names a file outside the scaffold also has: a package's own document, not a citation of the scaffold. */
const own = new Set(files.filter((f) => !f.startsWith(".seq/")).map((f) => basename(f)));
const cited = [...scaffold].filter((n) => !own.has(n));
const escape = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
const NAMED = new RegExp(`(?<![\\w.-])(${cited.map(escape).join("|")})(?![\\w-])`);
const PATH = /(?<![\w.-])\.seq\//;

const problems: string[] = [];
const dirty = new Set<string>();
for (const f of files) {
  if (!(f.startsWith("packages/") || f.startsWith("cli/")) || f === POINTER || !TEXT.test(f)) continue;
  const lines = readFileSync(join(ROOT, f), "utf8").split("\n");
  lines.forEach((line, i) => {
    const m = NAMED.exec(line) ?? PATH.exec(line);
    if (m && PENDING.has(f)) dirty.add(f);
    else if (m) problems.push(`${f}:${i + 1}: names the scaffold ("${m[0]}"): cite the spec section that holds the decision`);
  });
}
for (const f of PENDING) if (!dirty.has(f)) problems.push(`${f}: passes now: drop it from PENDING`);
finish("no-design-citations", problems);
