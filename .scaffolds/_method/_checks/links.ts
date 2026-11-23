// The link check: every relative Markdown link in a tracked file under `.scaffolds/` resolves to a tracked or present
// file or folder, and every `node .scaffolds/<path>.ts` command names a check that exists. A link inside a closed work
// folder (an `execution/<id>/` whose `state.md` ends with a `Closed` line) or a record (`records/`) is history, never
// edited after its close: reported as a count with `--history` listing them, not failed. A link or path with a
// `<placeholder>` in it (the method's templates) is skipped, as are links in code spans and fenced blocks, external
// links (`scheme:`), absolute paths and bare fragments; a fragment is not checked (the package's id check checks the
// package's).
// The method's: every scaffold is linked the same way, and one run covers all of `.scaffolds/`, so it takes no
// instance (`--instance` is accepted and narrows it to that scaffold, the package and the method).
// Usage: links.ts [root] [--instance .scaffolds/<name>] [--history]
import { existsSync, readFileSync } from "node:fs";
import { dirname, join, normalize } from "node:path";
import { args, finish, trackedFiles } from "../../principles/_checks/lib.ts";

const argv = process.argv.slice(2);
const HISTORY = argv.includes("--history");
process.argv.splice(2, process.argv.length, ...argv.filter((a) => a !== "--history"));
const { root: ROOT, instance: INSTANCE } = args();
const SCOPE = (f: string): boolean =>
  f.startsWith(".scaffolds/") &&
  (INSTANCE === undefined ||
    f.startsWith(`${INSTANCE}/`) ||
    !/^\.scaffolds\/[^_/][^/]*\//.test(f) ||
    /^\.scaffolds\/principles\//.test(f));

const tracked = trackedFiles(ROOT);
const present = new Set(tracked);
const dirs = new Set<string>();
for (const f of tracked) for (let d = dirname(f); d !== "." && !dirs.has(d); d = dirname(d)) dirs.add(d);
const exists = (p: string): boolean => present.has(p) || dirs.has(p) || existsSync(join(ROOT, p));

/** Closed work folders: `<scaffold>/execution/<id>/` whose `state.md` ends with a `Closed` line. */
const closed = new Set<string>();
for (const f of tracked) {
  const m = /^(\.scaffolds\/[^/]+\/execution\/[^/]+)\/state\.md$/.exec(f);
  if (!m || !existsSync(join(ROOT, f))) continue;
  const lines = readFileSync(join(ROOT, f), "utf8").split("\n").map((l) => l.trim()).filter(Boolean);
  if (/^Closed\b/.test(lines[lines.length - 1] ?? "")) closed.add(`${m[1]}/`);
}
const isHistory = (f: string): boolean =>
  /^\.scaffolds\/[^/]+\/records\//.test(f) || [...closed].some((c) => f.startsWith(c));

/** Relative link targets of a Markdown text, code spans and fenced blocks blanked first. */
function links(text: string): { target: string; line: number }[] {
  const lines = text.split("\n");
  let fence = false;
  const out: { target: string; line: number }[] = [];
  lines.forEach((raw, i) => {
    if (/^\s*(```|~~~)/.test(raw)) {
      fence = !fence;
      return;
    }
    if (fence) return;
    const line = raw.replace(/`[^`]*`/g, "``");
    for (const m of line.matchAll(/\]\(<([^>]+)>(?:\s+"[^"]*")?\)|\]\(([^)\s]+)(?:\s+"[^"]*")?\)|^\s*\[[^\]]+\]:\s*<?([^\s>]+)>?/g))
      out.push({ target: (m[1] ?? m[2] ?? m[3]) as string, line: i + 1 });
  });
  return out;
}

const problems: string[] = [];
const history: string[] = [];
let checked = 0;
for (const f of tracked) {
  if (!SCOPE(f) || !existsSync(join(ROOT, f))) continue;
  const isMd = f.endsWith(".md");
  if (!isMd && !/\.(ts|json)$/.test(f)) continue;
  const text = readFileSync(join(ROOT, f), "utf8");
  const found: string[] = [];
  if (isMd)
    for (const { target, line } of links(text)) {
      const path = target.split("#")[0] as string;
      if (!path || /^[a-z][a-z0-9+.-]*:/i.test(path) || path.startsWith("/") || /[<>]/.test(path)) continue;
      let decoded = path;
      try {
        decoded = decodeURIComponent(path);
      } catch {}
      checked++;
      const to = normalize(join(dirname(f), decoded)).replace(/\/$/, "");
      if (!exists(to)) found.push(`${f}:${line}: link to ${target} does not resolve (${to})`);
    }
  text.split("\n").forEach((line, i) => {
    for (const m of line.matchAll(/\bnode\s+(\.scaffolds\/[^\s`'"()]+\.ts)\b/g)) {
      const path = m[1] as string;
      if (/[<>]/.test(path)) continue;
      checked++;
      if (!exists(path)) found.push(`${f}:${i + 1}: command names ${path}, which does not exist`);
    }
  });
  (isHistory(f) ? history : problems).push(...found);
}
console.log(`links: ${checked} link(s) and command(s) checked; ${history.length} in closed history (reported, not failed)`);
if (HISTORY && history.length) console.log(history.map((h) => `  - ${h}`).join("\n"));
finish("links", problems);
