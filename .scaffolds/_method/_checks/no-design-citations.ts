// The no-design-citations check (§execution-unreferenced, carried to the whole scaffold: a scaffold is never cited by
// code, comments, tests or specs): nothing in the repository's code names a design or a source file of the scaffold
// (`<instance>/design/*.md`, `<instance>/sources/**/*.md`, read at run time) or a path into `.scaffolds/`; a decision
// is cited by the spec section that holds it. A file name passes only in a file of a package that holds a file of that
// name itself (its own `_spec/` documents), never because some other package has one.
// The method's, not an instance's: every scaffold is kept and never cited. The instance's part is
// `<instance>/_checks/no-design-citations.json`, required:
//   { "code": "<regex: the paths checked, e.g. the packages>", "own": "<regex: the scaffold's own packages>",
//     "pointers": ["<the files allowed to point at the scaffold>"], "pending": ["<files not yet clean>"] }
// The scaffold's own packages are held to their own documents; the other packages' documents are not the scaffold's to
// cite or to police: a name any of them holds passes there. `pending` lists files that a lane may not edit (the
// orchestrator's hot files) and only shrinks; a file leaves it when it is clean (the check then says so).
// Usage: no-design-citations.ts [root] --instance .scaffolds/<name> [--self-check] (self-check: a planted design file
// name in a comment of the scaffold's own package must fail, one in a package that holds that name itself must pass).
import { existsSync, readFileSync } from "node:fs";
import { basename, join } from "node:path";
import { finish, instanceFile, methodArgs, trackedFiles, usage } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, flags } = methodArgs({ flags: ["--self-check"] });
const SELF = flags.has("--self-check");
const config = instanceFile<{ code?: string; own?: string; pointers?: string[]; pending?: string[] }>(
  ROOT,
  INSTANCE,
  "no-design-citations.json",
);
if (!config.code || !config.own) usage(`${INSTANCE}/_checks/no-design-citations.json: "code" and "own" are required`);
const CODE = new RegExp(config.code);
const OWN = new RegExp(config.own);
const POINTERS = new Set(config.pointers ?? []);
const PENDING = new Set(config.pending ?? []);
const TEXT = /\.(ts|tsx|mts|mjs|js|json|md|mdx|css|scss|html|yml|yaml|txt)$/;
const SCAFFOLD = `${INSTANCE}/`;

const files = trackedFiles(ROOT).filter((f) => existsSync(join(ROOT, f)));

/** The scaffold's design and source file names (without the folder, README.md left out: every folder has one). */
const scaffold = new Set<string>();
for (const f of files) {
  if (!f.endsWith(".md") || basename(f) === "README.md") continue;
  if (f.startsWith(`${SCAFFOLD}design/`) || f.startsWith(`${SCAFFOLD}sources/`)) scaffold.add(basename(f));
}
/** The package folder of a file: the nearest ancestor holding a tracked `package.json`, else its first two folders. */
const roots = new Set(files.filter((f) => f.endsWith("/package.json")).map((f) => f.slice(0, -"/package.json".length)));
const packageOf = (f: string): string => {
  for (let d = f.slice(0, Math.max(0, f.lastIndexOf("/"))); d !== ""; d = d.slice(0, Math.max(0, d.lastIndexOf("/"))))
    if (roots.has(d)) return d;
  return f.split("/").slice(0, 2).join("/");
};
/** The names each package holds outside the scaffolds: its own documents, not citations of the scaffold. */
const held = new Map<string, Set<string>>();
for (const f of files) {
  if (f.startsWith(".scaffolds/")) continue;
  const p = packageOf(f);
  if (!held.has(p)) held.set(p, new Set());
  held.get(p)?.add(basename(f));
}
const escape = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
const NAMED = scaffold.size
  ? new RegExp(`(?<![\\w.-])(${[...scaffold].map(escape).join("|")})(?![\\w-])`, "g")
  : undefined;
const PATH = /(?<![\w.-])\.scaffolds\//;

const outside = new Set(files.filter((f) => !f.startsWith(".scaffolds/") && !OWN.test(f)).map((f) => basename(f)));

/** What a line of the file `f` names of the scaffold, if anything. */
function names(f: string, line: string): string | undefined {
  const own = OWN.test(f) ? held.get(packageOf(f)) : outside;
  if (NAMED) for (const m of line.matchAll(NAMED)) if (!own?.has(m[1] as string)) return m[0];
  return PATH.exec(line)?.[0];
}

if (SELF) {
  const problems: string[] = [];
  // a name some other package holds too tests both halves; any name tests the first
  const otherHeld = (n: string) => [...held].some(([p, h]) => h.has(n) && !OWN.test(`${p}/x`));
  const name = [...scaffold].sort().find(otherHeld) ?? [...scaffold].sort()[0];
  const sample = files.find((f) => OWN.test(f) && CODE.test(f) && TEXT.test(f));
  if (name === undefined || sample === undefined) usage("--self-check needs a design or source file and an own package");
  const planted = `// the rule (${name}, Consequences)`;
  // a planted file in the scaffold's own package, unless that package holds the name itself
  const ownPackage = packageOf(sample);
  if (!held.get(ownPackage)?.has(name) && names(`${ownPackage}/src/x.ts`, planted) === undefined)
    problems.push(`a planted ${name} in ${ownPackage} passed`);
  const holder = [...held].find(([p, n]) => n.has(name) && !OWN.test(`${p}/x`));
  if (holder && names(`${holder[0]}/src/x.ts`, planted) !== undefined)
    problems.push(`${holder[0]} holds a ${name} of its own and the citation failed there`);
  finish("no-design-citations --self-check", problems);
}

const problems: string[] = [];
const dirty = new Set<string>();
for (const f of files) {
  if (!CODE.test(f) || POINTERS.has(f) || !TEXT.test(f)) continue;
  const lines = readFileSync(join(ROOT, f), "utf8").split("\n");
  lines.forEach((line, i) => {
    const m = names(f, line);
    if (m && PENDING.has(f)) dirty.add(f);
    else if (m) problems.push(`${f}:${i + 1}: names the scaffold ("${m}"): cite the spec section that holds the decision`);
  });
}
for (const f of PENDING) if (!dirty.has(f)) problems.push(`${f}: passes now: drop it from no-design-citations.json's pending`);
finish("no-design-citations", problems);
