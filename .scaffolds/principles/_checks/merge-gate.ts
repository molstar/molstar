// The merge check (§execution-unreferenced; README.md here): no tracked file outside the work folder refers to a path
// inside it. The folder itself is committed and kept as the work's archive (§permanent-record), and its files may refer
// to each other. The work folder is the instance's `execution/` when the project has an instance, else `.execution/` at
// the root (the package's role name); any `.execution/` folder counts.
// A reference is a root-relative path in any text, or a relative Markdown link that resolves into the folder.
// Rule and docs text that talks about the folder in general is allowlisted: the package, the instance's principles
// and skills, `.agents/`, every AGENTS.md, plus `--allow <regex>` for the project's own.
// Usage: merge-gate.ts [root] [--package .principles] [--instance .<name>] [--allow <regex>]...
import { readFileSync } from "node:fs";
import { dirname, join, normalize } from "node:path";
import { args, finish, inWorkFolder, isLockfile, trackedFiles, workFolder } from "./lib.ts";

const { root: ROOT, pkg: PKG, instance: INSTANCE, allow } = args({ allow: true });
const esc = (s: string): string => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
const WORK = workFolder(INSTANCE);
const ALLOWED = [
  new RegExp(`^${esc(PKG)}/`),
  ...(INSTANCE ? [new RegExp(`^${esc(INSTANCE)}/(?:principles|_skills)/`)] : []),
  /^\.agents\//,
  /(^|\/)AGENTS\.md$/,
  ...allow,
];
/** A path inside the folder (not the bare word or the folder itself), not the tail of a longer name. */
const PATH_REFERENCE = new RegExp(`(?<![\\w.-])(?:${esc(WORK)}|\\.execution/)[\\w.-]`);
const inWork = (p: string): boolean => inWorkFolder(p, INSTANCE);
/** Relative Markdown links (inline and reference-style) that land inside the work folder. */
const linksIntoWork = (f: string, text: string): string[] =>
  [
    ...text.matchAll(
      /\]\(<([^>\n]+)>(?:\s+"[^"]*")?\)|\]\(([^)\s<]+)(?:\s+"[^"]*")?\)|^\s*\[[^\]]+\]:\s*<?([^\s>]+)>?/gm,
    ),
  ]
    .map((m) => (m[1] ?? m[2] ?? m[3] ?? "").split("#")[0] as string)
    .filter((l) => l && !/^[a-z][a-z0-9+.-]*:/i.test(l) && !l.startsWith("/"))
    .filter((l) => inWork(normalize(join(dirname(f), decodeURIComponent(l)))));

const problems: string[] = [];
for (const f of trackedFiles(ROOT)) {
  if (inWork(f)) continue; // the archive refers to itself
  if (ALLOWED.some((re) => re.test(f)) || isLockfile(f)) continue;
  let text: string;
  try {
    text = readFileSync(join(ROOT, f), "utf8");
  } catch {
    continue; // deleted in the working tree, or a link to a folder
  }
  if (text.includes("\0")) continue; // binary
  if (PATH_REFERENCE.test(text))
    problems.push(`${f} refers to a path inside the work folder (§execution-unreferenced)`);
  else if (/\.mdx?$/.test(f))
    for (const l of linksIntoWork(f, text))
      problems.push(`${f} links into the work folder: ${l} (§execution-unreferenced)`);
}

finish("merge-gate", problems);
