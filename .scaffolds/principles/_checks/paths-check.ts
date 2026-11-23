// The path check (§machine-invariant; README.md here): no tracked file names a path on one person's machine (a home
// folder, a temp folder, a drive path). References point inside the repository (relative paths), at a published
// location (a repository URL), or name the thing in words. The work folder is committed and checked like the rest.
// Caught: a home folder (macOS Users/<name>, Linux home/<name>/...), a temp folder (tmp, private/tmp, var/folders,
// at the filesystem root) or a drive path (a letter, a colon, a slash or backslash), also in a file URL or
// JSON-escaped.
// Usage: paths-check.ts [root] [--allow <regex>]...
import { readFileSync } from "node:fs";
import { join } from "node:path";
import { args, finish, isLockfile, trackedFiles } from "./lib.ts";

const { root: ROOT, allow } = args({ allow: true });
/** After a start, a space, punctuation or `file://` (not inside a URL path); a Linux home needs a segment after the
 * name, so an app route like `/home/settings` is not one. */
const MACHINE_PATH =
  /(^|[\s`"'(=:[,<>|;{]|file:\/\/)(\/Users\/[^/\s`"')<>|;]+|\/home\/[^/\s`"')<>|;]+(?=\/)|\/(?:private\/)?tmp\/|\/private\/var\/|\/var\/folders\/|[A-Za-z]:(?:\\{1,2}|\/)(?=[\w.~-]))/m;

const files = trackedFiles(ROOT);
const problems: string[] = [];
for (const f of files) {
  if (isLockfile(f) || allow.some((re) => re.test(f))) continue;
  let text: string;
  try {
    text = readFileSync(join(ROOT, f), "utf8");
  } catch {
    continue; // deleted in the working tree
  }
  if (text.includes("\0")) continue; // binary
  for (const m of text.matchAll(new RegExp(MACHINE_PATH.source, "gm"))) {
    const line = text.slice(0, m.index).split("\n").length;
    problems.push(`${f}:${line} names a path on one machine (${m[2]}…); use a repository-relative path or a URL`);
  }
}
finish("paths-check", problems);
