// Shared by the package's checks (README.md here). Node only, no dependencies.
import { execFileSync } from "node:child_process";
import { existsSync, statSync } from "node:fs";
import { isAbsolute, join, normalize, relative, resolve } from "node:path";

export interface CheckArgs {
  /** The repository root, absolute; the one positional argument, else the working directory. */
  root: string;
  /** The package folder, relative to the root (`--package`, default `.principles`). */
  pkg: string;
  /** The project's instance folder, relative to the root (`--instance`), if it has one. */
  instance: string | undefined;
  /** Extra allowlist patterns over file paths (`--allow <regex>`, repeatable), for checks that take one. */
  allow: RegExp[];
}

/** What a check accepts: `allow`, the `--allow` flag; `needsPackage`, the package folder must exist even by default. */
export interface CheckOptions {
  allow?: boolean;
  needsPackage?: boolean;
}

/** Usage of every check: `<check>.ts [root] [--package <dir>] [--instance <dir>] [--allow <regex>]...` */
export function parseArgs(argv: readonly string[] = process.argv.slice(2), opts: CheckOptions = {}): CheckArgs {
  let rootArg: string | undefined;
  let pkgArg: string | undefined;
  let instanceArg: string | undefined;
  const allow: RegExp[] = [];
  const value = (flag: string, v: string | undefined): string => {
    if (v === undefined || v === "" || v.startsWith("--")) throw new Error(`${flag} needs a value`);
    return v;
  };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i] as string;
    if (a === "--package") pkgArg = value(a, argv[++i]);
    else if (a === "--instance") instanceArg = value(a, argv[++i]);
    else if (a === "--allow") {
      if (!opts.allow) throw new Error("--allow is not an option of this check");
      allow.push(new RegExp(value(a, argv[++i])));
    } else if (a.startsWith("-")) throw new Error(`unknown flag ${a}`);
    else if (rootArg !== undefined) throw new Error(`one root only (got ${rootArg} and ${a})`);
    else rootArg = a;
  }
  const root = resolve(rootArg ?? ".");
  if (!existsSync(root) || !statSync(root).isDirectory()) throw new Error(`root ${root} is not a folder`);
  /** A folder inside the root, as a clean relative path; an absolute one is made relative, one outside is refused. */
  const dir = (flag: string, v: string, mustExist = true): string => {
    const rel = normalize(relative(root, resolve(root, v))).replace(/\/$/, "");
    if (rel === "" || rel === "." || rel === ".." || rel.startsWith("../") || isAbsolute(rel))
      throw new Error(`${flag} ${v} must be a folder inside the root`);
    if (mustExist && (!existsSync(join(root, rel)) || !statSync(join(root, rel)).isDirectory()))
      throw new Error(`${flag} ${v}: no such folder under ${root}`);
    return rel;
  };
  return {
    root,
    pkg: dir("--package", pkgArg ?? ".principles", pkgArg !== undefined || opts.needsPackage === true),
    instance: instanceArg === undefined ? undefined : dir("--instance", instanceArg),
    allow,
  };
}

/** Runs a check's arguments through parseArgs; a usage error exits 2 with the message (never a silent pass). */
export function args(opts: CheckOptions = {}): CheckArgs {
  try {
    return parseArgs(process.argv.slice(2), opts);
  } catch (e) {
    console.error(`usage: ${(e as Error).message}`);
    process.exit(2);
  }
}

export const git = (cwd: string, ...a: string[]): string =>
  execFileSync("git", a, { cwd, encoding: "utf8", maxBuffer: 64 * 1024 * 1024 });

/** Tracked files, NUL-separated so names with non-ASCII or special characters come back as written, not quoted. */
export const trackedFiles = (root: string): string[] => git(root, "ls-files", "-z").split("\0").filter(Boolean);

/** The work folder (§execution-unreferenced): the instance's `execution/`, else `.execution/` at the root; any
 * `.execution/` folder counts. Committed and kept as each work's archive (§permanent-record). */
export const workFolder = (instance: string | undefined): string =>
  instance ? `${instance}/execution/` : ".execution/";
export const inWorkFolder = (f: string, instance: string | undefined): boolean =>
  f.startsWith(workFolder(instance)) || f.startsWith(".execution/") || f.includes("/.execution/");

/** Lockfiles are generated, large and full of hashes: no check reads them. */
export const isLockfile = (f: string): boolean =>
  /(^|\/)(pnpm-lock\.yaml|package-lock\.json|yarn\.lock|bun\.lockb?|Cargo\.lock|uv\.lock|poetry\.lock|go\.sum)$/.test(
    f,
  );

/** Prints the problems and exits 1 when there are any; every check's last step. */
export function finish(name: string, problems: readonly string[]): never {
  if (problems.length) {
    console.error(`${name}: ${problems.length} problem(s)\n${problems.map((p) => `  - ${p}`).join("\n")}`);
    process.exit(1);
  }
  console.log(`${name}: ok`);
  process.exit(0);
}
