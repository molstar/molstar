// The id check (_core/system.md, "Ids" and "Code patterns"; README.md here). Every rule item of a rule set and every
// pattern README carries a `{#id}` anchor. A code pattern lives in <module>/bootstrap-patterns/<id>/ in the package, or
// <module>/stable-patterns/<id>/ or sketch-patterns/<id>/ in the instance (status is the folder), named by its id.
// Ids are unique across the package and the instance. Every citation in a tracked file resolves (PLACEHOLDERS
// excepted): the section sign followed by a kebab-case id. An old number (letter, section, optional item) is rejected
// with a hint to cite the id; anything else after the sign is malformed. Tombstones that only redirect are reported
// (not failed), so their citations get updated. Relative links (pages, images, #fragments) resolve in the package, the
// instance's principles, design, sources and skills, AGENTS.md and CLAUDE.md. The package is self-contained: a link
// from it never leaves it, and a citation in it resolves inside it; _core cites only _core.
// The work folder is not checked: it is an archive, never cited and never edited after the close (§permanent-record).
// `--allow <regex>` skips host files whose use of the section sign is not a citation.
// Usage: principles-ids.ts [root] [--package .principles] [--instance .<name>] [--allow <regex>]...
import { existsSync, readFileSync, statSync } from "node:fs";
import { dirname, join, relative } from "node:path";
import { args, finish, inWorkFolder, isLockfile, trackedFiles } from "./lib.ts";

const { root: ROOT, pkg: PKG, instance: INSTANCE, allow } = args({ allow: true, needsPackage: true });
const TRACKED = trackedFiles(ROOT);
const BASE = join(ROOT, PKG);
const IN_PKG = (rel: string): boolean => rel.startsWith(`${PKG}/`);
/** The package and the instance: modules in both define ids; rule sets are laid out the same; pattern folders differ (below). */
const BASES = [BASE, ...(INSTANCE ? [join(ROOT, INSTANCE, "principles")] : [])];
/** A package rule's word budget (§rule-budget): a rule states what and why; the how and the evidence live elsewhere. */
const RULE_WORDS = 180;
const S = "§"; // the section sign; no citation in this file is followed by an id
/** _core stands alone (system.md, "Adapting"): what it and its skill cite must be defined in _core. */
const STANDALONE_DIRS = ["_core/", "_skills/principles/"].map((d) => `${PKG}/${d}`);
const standalone = (rel: string): boolean => STANDALONE_DIRS.some((d) => rel.startsWith(d));

/** Known false positives: illustrative links inside examples, by file name. */
const ALLOWED_LINKS: Record<string, readonly string[]> = { "spec-mirror": ["module.md", "module.md#documentstore"] };

/** Tracked Markdown files under a folder (relative to the root), as absolute paths: what CI sees, not local drafts. */
const trackedMd = (dir: string): string[] =>
  TRACKED.filter((f) => f.startsWith(`${dir}/`) && f.endsWith(".md")).map((f) => join(ROOT, f));
/** True when `target` is outside `base` (a sibling named `..x` is inside). */
const leaves = (base: string, target: string): boolean => {
  const r = relative(base, target);
  return r === ".." || r.startsWith("../") || r.startsWith("/");
};

const ID = "[a-z][a-z0-9]*(?:-[a-z0-9]+)*";
const ANCHOR = new RegExp(`\\{#(${ID})\\}`, "g");
const HAS_ANCHOR = new RegExp(`\\{#(${ID})\\}`);
const CITE = new RegExp(`${S}([A-Za-z0-9][\\w.-]*)`, "g");
const NAME = new RegExp(`^(${ID})`);
const NUMBER = /^[A-Z]\d+(?:\.\d+)?(?![\w])/;
/** A code pattern's folder is its status: bootstrap in the package; stable or sketch in the instance. */
const PACKAGE_PATTERN_DIRS = new Set(["bootstrap-patterns"]);
const INSTANCE_PATTERN_DIRS = new Set(["stable-patterns", "sketch-patterns"]);
const LANGUAGE_FILE = /^(typescript|javascript|python|rust|go|java|kotlin|swift|csharp|cpp|c)\.md$/;
/** Placeholder ids that docs use to show the citation form; they are examples, not citations. */
const PLACEHOLDERS = new Set(["id", "kebab-id", "old-id", "new-id"]);
const REDIRECT = new RegExp(`\\{#(${ID})\\}\\s*\\*(?:moved to|see) ${S}(${ID})\\*\\s*$`);

const problems: string[] = [];
const notices: string[] = [];
const defs = new Map<string, string>(); // id -> file that defines it
const redirects = new Map<string, string>(); // tombstone id -> target id

/** Lines outside fenced code, with inline code removed (anchors, items and links there are examples). */
const proseLines = (text: string): { line: string; raw: string; n: number }[] => {
  const out: { line: string; raw: string; n: number }[] = [];
  let fenced = false;
  text.split("\n").forEach((raw, i) => {
    if (/^\s*```/.test(raw)) fenced = !fenced;
    else if (!fenced) out.push({ line: raw.replace(/`[^`]*`/g, ""), raw, n: i + 1 });
  });
  return out;
};

const principleFiles = BASES.flatMap((b) => trackedMd(relative(ROOT, b)).map((f) => ({ f, base: b })));
const patternDirs = new Map<string, boolean>(); // <base>/<module>/<pattern-dir>/<id> -> has a README
for (const { f, base } of principleFiles) {
  const rel = relative(ROOT, f);
  const parts = relative(base, f).split("/");
  const isRuleSet =
    (parts.length === 2 && parts[1] === "rules.md") ||
    rel === `${PKG}/_core/workflow.md` ||
    rel === `${PKG}/_core/lanes.md`;
  const budgeted = isRuleSet && base === BASE; // the package's rules have a word budget; an instance sets its own
  let item: { id: string; n: number; words: number } | undefined;
  const closeItem = (): void => {
    if (item && item.words > RULE_WORDS)
      problems.push(
        `${rel}:${item.n}: rule {#${item.id}} has ${item.words} words, over ${RULE_WORDS} (§rule-budget: move the how to a pattern, the evidence to a record)`,
      );
    item = undefined;
  };
  const allowed = base === BASE ? PACKAGE_PATTERN_DIRS : INSTANCE_PATTERN_DIRS;
  const other = base === BASE ? INSTANCE_PATTERN_DIRS : PACKAGE_PATTERN_DIRS;
  const inPatterns = parts[0] !== "_template" && allowed.has(parts[1] ?? "");
  if (parts[0] !== "_template" && parts.length > 2 && other.has(parts[1] ?? ""))
    problems.push(
      `${rel}: ${parts[1]}/ belongs ${base === BASE ? "in an instance" : "in the package"}; here patterns go in ${[...allowed].join(" or ")}/`,
    );
  const isPattern = inPatterns && parts.length === 4 && parts[3] === "README.md";
  if (inPatterns && parts.length >= 4) {
    const dir = [relative(ROOT, base), ...parts.slice(0, 3)].join("/");
    patternDirs.set(dir, (patternDirs.get(dir) ?? false) || isPattern);
  }
  if (!existsSync(f)) continue; // tracked, deleted in the working tree
  if (parts[0] !== "_template" && LANGUAGE_FILE.test(parts.at(-1) ?? "") && !(inPatterns && parts.length === 4))
    problems.push(`${rel}: an implementation file outside a pattern folder (<module>/${[...allowed].join("|")}/<id>/)`);
  let section = "";
  for (const { line, raw, n } of proseLines(readFileSync(f, "utf8"))) {
    for (const m of line.matchAll(ANCHOR)) {
      const id = m[1] as string;
      const prev = defs.get(id);
      if (prev) problems.push(`${rel}: id {#${id}} is already defined in ${prev}`);
      else defs.set(id, rel);
    }
    const r = REDIRECT.exec(line.replace(/^\s*(?:[-*]|\d+\.)\s+/, ""));
    if (r) redirects.set(r[1] as string, r[2] as string);
    if (budgeted) {
      if (/^\d+\. /.test(raw) || raw.startsWith("#") || !raw.trim()) closeItem();
      const anchor = /^\d+\. /.test(raw) && section !== "## In brief" ? HAS_ANCHOR.exec(line)?.[1] : undefined;
      if (anchor) item = { id: anchor, n, words: 0 };
      if (item) item.words += raw.split(/\s+/).filter(Boolean).length;
    }
    if (isRuleSet) {
      if (raw.startsWith("## ")) {
        section = raw.trim();
        if (section !== "## In brief" && !HAS_ANCHOR.test(line)) problems.push(`${rel}:${n}: section has no {#id}`);
      } else if (/^\d+\. /.test(raw) && section && section !== "## In brief" && !HAS_ANCHOR.test(line))
        problems.push(`${rel}:${n}: rule has no {#id} anchor`);
    }
    if (isPattern && n === 1) {
      const id = HAS_ANCHOR.exec(line)?.[1];
      if (!id) problems.push(`${rel}: pattern title has no {#id} anchor`);
      else if (id !== parts[2]) problems.push(`${rel}: pattern id {#${id}} differs from its folder name ${parts[2]}`);
    }
  }
  closeItem();
}

for (const [dir, hasReadme] of patternDirs)
  if (!hasReadme) problems.push(`${dir}/: a pattern folder without a README.md (its id, concept and cases)`);

// GitHub-style heading slugs (a repeated heading gets -1, -2, ...), plus explicit anchors, per Markdown file
const fragmentsOf = new Map<string, Set<string>>();
const fragments = (file: string): Set<string> => {
  let set = fragmentsOf.get(file);
  if (!set) {
    set = new Set();
    const slugCount = new Map<string, number>();
    for (const { line, raw } of proseLines(readFileSync(file, "utf8"))) {
      for (const m of line.matchAll(ANCHOR)) set.add(m[1] as string);
      const h = /^#{1,6} (.*)$/.exec(raw);
      if (h) {
        const slug = (h[1] as string)
          .replace(/\s*\{#[^}]+\}/g, "")
          .trim()
          .toLowerCase()
          .replace(/[^\p{L}\p{N}_\- ]/gu, "")
          .replace(/ /g, "-");
        const seen = slugCount.get(slug) ?? 0;
        slugCount.set(slug, seen + 1);
        set.add(seen ? `${slug}-${seen}` : slug);
      }
    }
    fragmentsOf.set(file, set);
  }
  return set;
};

const linkFiles = [
  ...principleFiles.map((p) => p.f),
  ...(INSTANCE ? ["design", "sources", "_skills"].flatMap((d) => trackedMd(`${INSTANCE}/${d}`)) : []),
  ...["AGENTS.md", "CLAUDE.md"].filter((f) => TRACKED.includes(f)).map((f) => join(ROOT, f)),
].filter((f) => existsSync(f));
/** Inline links (`](x)`, `](<x y>)`, `](x "title")`) and reference definitions (`[r]: x`). */
const LINK = /\]\(<([^>\n]+)>(?:\s+"[^"]*")?\)|\]\(([^)\s<]+)(?:\s+"[^"]*")?\)|^\s*\[[^\]]+\]:\s*<?([^\s>]+)>?/g;
for (const f of linkFiles) {
  const rel = relative(ROOT, f);
  for (const { line, n } of proseLines(readFileSync(f, "utf8")))
    for (const m of line.matchAll(LINK)) {
      const link = (m[1] ?? m[2] ?? m[3]) as string;
      if (/^[a-z][a-z0-9+.-]*:/i.test(link)) continue;
      if (Object.entries(ALLOWED_LINKS).some(([k, v]) => rel.includes(k) && v.includes(link))) continue;
      const [path = "", frag] = link.split("#");
      const target = !path ? f : join(path.startsWith("/") ? ROOT : dirname(f), decodeURIComponent(path));
      if (IN_PKG(rel) && leaves(BASE, target))
        problems.push(`${rel}:${n}: link ${link} leaves ${PKG}/ (the package is self-contained)`);
      else if (!existsSync(target)) problems.push(`${rel}:${n}: broken link ${link}`);
      else if (frag && target.endsWith(".md") && statSync(target).isFile() && !fragments(target).has(frag))
        problems.push(`${rel}:${n}: broken fragment ${link}`);
    }
}

const tracked = TRACKED;
const cited = new Map<string, number>();
for (const f of tracked) {
  // the work folder is an archive: its citations name ids as they were then, and it is never edited after the close
  if (isLockfile(f) || inWorkFolder(f, INSTANCE) || allow.some((re) => re.test(f))) continue;
  let text: string;
  try {
    text = readFileSync(join(ROOT, f), "utf8");
  } catch {
    continue; // deleted in the working tree, or a link to a folder
  }
  if (text.includes("\0") || !text.includes(S)) continue;
  const lineOf = (i: number) => text.slice(0, i).split("\n").length;
  for (const m of text.matchAll(CITE)) {
    const token = m[1] as string;
    const where = `${f}:${lineOf(m.index ?? 0)}`;
    const name = NAME.exec(token)?.[1];
    const num = NUMBER.exec(token)?.[0];
    if (name && PLACEHOLDERS.has(name)) continue;
    if (name) {
      cited.set(name, (cited.get(name) ?? 0) + 1);
      if (!defs.has(name)) problems.push(`${where}: unresolved citation ${S}${name}`);
      else if (standalone(f) && !standalone(defs.get(name) as string))
        problems.push(`${where}: ${S}${name} is defined outside _core; _core cites no other module`);
      else if (IN_PKG(f) && !IN_PKG(defs.get(name) as string))
        problems.push(`${where}: ${S}${name} is defined in the instance; the package cites only itself`);
    } else if (num) {
      problems.push(`${where}: old numbered citation ${S}${num}; cite the rule by its id`);
    } else problems.push(`${where}: malformed citation ${S}${token} (an id is kebab-case)`);
  }
}

for (const [id, target] of redirects) {
  if (!defs.has(target)) problems.push(`tombstone ${S}${id} redirects to unknown ${S}${target}`);
  notices.push(`tombstone ${S}${id} only redirects to ${S}${target} (${cited.get(id) ?? 0} citation(s) to update)`);
}
if (notices.length)
  console.log(`principles-ids: ${notices.length} notice(s)\n${notices.map((n) => `  - ${n}`).join("\n")}`);
console.log(`principles-ids: ${defs.size} ids`);
finish("principles-ids", problems);
