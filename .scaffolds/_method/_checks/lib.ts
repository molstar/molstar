// Shared by the method's checks (../README.md, Checks): the principles package's arguments with `--instance`
// required, each check's own flags, and the instance's file a check reads its scaffold-specific part from
// (`<instance>/_checks/<file>`). Node only, no dependencies. Cites only the package's ids: the method serves every
// scaffold, so an instance's id would not resolve in another.
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { type CheckArgs, args } from "../../principles/_checks/lib.ts";

export { finish, trackedFiles } from "../../principles/_checks/lib.ts";

export interface MethodArgs extends CheckArgs {
  /** The scaffold checked, relative to the root (`--instance .scaffolds/<name>`, required). */
  instance: string;
  /** The check's boolean flags that were given. */
  flags: Set<string>;
  /** The check's valued options that were given. */
  values: Map<string, string>;
}

export function usage(message: string): never {
  console.error(`usage: ${message}`);
  process.exit(2);
}

/** `<check>.ts [root] --instance <dir> [flags]`: a bad or missing argument exits 2, never a silent pass. */
export function methodArgs(own: { flags?: string[]; values?: string[] } = {}): MethodArgs {
  const argv = process.argv.slice(2);
  const flags = new Set<string>();
  const values = new Map<string, string>();
  const rest: string[] = [];
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i] as string;
    if (own.flags?.includes(a)) flags.add(a);
    else if (own.values?.includes(a)) {
      const v = argv[++i];
      if (v === undefined || v === "" || v.startsWith("--")) usage(`${a} needs a value`);
      values.set(a, v);
    } else rest.push(a);
  }
  process.argv.splice(2, process.argv.length, ...rest);
  const parsed = args();
  if (parsed.instance === undefined) usage("--instance .scaffolds/<name> is required");
  return { ...parsed, instance: parsed.instance, flags, values };
}

/** The instance's file for a check, parsed as JSON: `fallback` when it is absent, else a usage error naming it. */
export function instanceFile<T>(root: string, instance: string, file: string, fallback?: T): T {
  const path = join(root, instance, "_checks", file);
  if (!existsSync(path)) {
    if (fallback !== undefined) return fallback;
    usage(`${instance}/_checks/${file} is missing: the check reads its scaffold's part from it`);
  }
  try {
    return JSON.parse(readFileSync(path, "utf8")) as T;
  } catch (e) {
    usage(`${instance}/_checks/${file}: ${(e as Error).message}`);
  }
}

/** The text of a Markdown file's `## <title>` section, or undefined when it has none. */
export function section(text: string, title: RegExp): string | undefined {
  const lines = text.split("\n");
  const from = lines.findIndex((l) => /^## /.test(l) && title.test(l.slice(3)));
  if (from < 0) return undefined;
  let to = lines.findIndex((l, i) => i > from && /^## /.test(l));
  if (to < 0) to = lines.length;
  return lines.slice(from + 1, to).join("\n");
}
