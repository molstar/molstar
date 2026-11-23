// The no-design-citations check (stage 1b of the sixth slice): this folder is temporary, so nothing under `packages/`
// (or the command-line packages in `cli/`) names a design or a source file of the scaffold (`.seq/design/*.md`,
// `.seq/sources/**/*.md`, read at run time) or a path into `.seq/`; a decision is cited by the spec section that holds
// it. A file name passes only in a file of a package that holds a file of that name itself (its own `_spec/`
// documents), never because some other package has one. The one pointer to the scaffold is
// `packages/seq/AGENTS.md`, until the scaffold is deleted.
// Usage: no-design-citations.ts [root] [--self-check] (self-check: a planted `sources.md` in a `packages/seq` comment
// must fail, one in a package that has its own `sources.md` must pass).
import { existsSync, readFileSync } from "node:fs";
import { basename, join } from "node:path";
import { args, finish, trackedFiles } from "../oz/_checks/lib.ts";

const SELF = process.argv.includes("--self-check");
if (SELF) process.argv.splice(process.argv.indexOf("--self-check"), 1);
const { root: ROOT } = args();
const POINTER = "packages/seq/AGENTS.md";
/** Files that name the scaffold when the check was written (2026-11-23) and that a lane may not edit (the
 * orchestrator's hot files): the list only shrinks, a file leaves it when it is clean (the check then says so). */
const PENDING = new Set<string>([]);
const TEXT = /\.(ts|tsx|mts|mjs|js|json|md|mdx|css|scss|html|yml|yaml|txt)$/;

const files = trackedFiles(ROOT).filter((f) => existsSync(join(ROOT, f)));

/** The scaffold's design and source file names (without the folder, README.md left out: every folder has one). */
const scaffold = new Set<string>();
for (const f of files) {
  if (!f.endsWith(".md") || basename(f) === "README.md") continue;
  if (f.startsWith(".seq/design/") || f.startsWith(".seq/sources/")) scaffold.add(basename(f));
}
/** The package folder of a file: the nearest ancestor holding a tracked `package.json`, else its first two folders. */
const roots = new Set(files.filter((f) => f.endsWith("/package.json")).map((f) => f.slice(0, -"/package.json".length)));
const packageOf = (f: string): string => {
  for (let d = f.slice(0, Math.max(0, f.lastIndexOf("/"))); d !== ""; d = d.slice(0, Math.max(0, d.lastIndexOf("/"))))
    if (roots.has(d)) return d;
  return f.split("/").slice(0, 2).join("/");
};
/** The names each package holds outside the scaffold: its own documents, not citations of the scaffold. */
const held = new Map<string, Set<string>>();
for (const f of files) {
  if (f.startsWith(".seq/")) continue;
  const p = packageOf(f);
  if (!held.has(p)) held.set(p, new Set());
  held.get(p)?.add(basename(f));
}
const escape = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
const NAMED = new RegExp(`(?<![\\w.-])(${[...scaffold].map(escape).join("|")})(?![\\w-])`, "g");
const PATH = /(?<![\w.-])\.seq\//;

/** Seq*'s own packages are held to their own documents; the other packages' documents (Prism's, the sequences) are
 * not Seq*'s to cite or to police: a name any of them holds passes there. */
const SEQ = /^(packages\/seq\/|cli\/seqx)/;
const outside = new Set(files.filter((f) => !f.startsWith(".seq/") && !SEQ.test(f)).map((f) => basename(f)));

/** What a line of the file `f` names of the scaffold, if anything. */
function names(f: string, line: string): string | undefined {
  const own = SEQ.test(f) ? held.get(packageOf(f)) : outside;
  for (const m of line.matchAll(NAMED)) if (!own?.has(m[1] as string)) return m[0];
  return PATH.exec(line)?.[0];
}

if (SELF) {
  const planted = "// the rule (sources.md, Consequences)";
  const bad = names("packages/seq/core/src/x.ts", planted);
  const own = [...held].find(([, n]) => n.has("sources.md"));
  const problems: string[] = [];
  if (bad === undefined && scaffold.has("sources.md")) problems.push("a planted sources.md in packages/seq passed");
  if (own && !own[0].startsWith("packages/seq") && names(`${own[0]}/src/x.ts`, planted) !== undefined)
    problems.push(`${own[0]} holds a sources.md of its own and the citation failed there`);
  finish("no-design-citations --self-check", problems);
  process.exit(0);
}

const problems: string[] = [];
const dirty = new Set<string>();
for (const f of files) {
  if (!(f.startsWith("packages/") || f.startsWith("cli/")) || f === POINTER || !TEXT.test(f)) continue;
  const lines = readFileSync(join(ROOT, f), "utf8").split("\n");
  lines.forEach((line, i) => {
    const m = names(f, line);
    if (m && PENDING.has(f)) dirty.add(f);
    else if (m) problems.push(`${f}:${i + 1}: names the scaffold ("${m}"): cite the spec section that holds the decision`);
  });
}
for (const f of PENDING) if (!dirty.has(f)) problems.push(`${f}: passes now: drop it from PENDING`);
finish("no-design-citations", problems);
