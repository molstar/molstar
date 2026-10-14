/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq* export checks (RB0).
//   The guard: no name is exported by two entries (of one package or of two) as different declarations, in
//   the type namespace or in the value namespace, because a host imports the packages together. Walked with the
//   TypeScript checker over every entry of every packages/sequences package; a name that resolves to the same
//   declaration (a re-export) is allowed, and so is a package's `./react` against its `./dom` (one slot contract).
//   Exit 1 when a name collides.
//   The contract check: every Seq* type a public signature names (an exported interface's, type's, class's public
//   members, a function's parameters and result) is exported by some entry, so a host or a renderer can name it. A
//   member of an exported namespace counts as exported; `/testing` entries are test helpers, not the contract. Checked
//   for every package, the core included; `--contract` lists the gaps without failing. Exit 1 on a gap.
//   The listing (--unused): the exported names no consumer mentions, for SW1 to drop or justify. A consumer is
//   any tracked .ts or .tsx file outside the entry's own package and the packages' test folders. A listing, not a
//   check: a name only the package's own tests use is listed too.
// Usage: node packages/sequences/scripts/exports.ts [--unused | --contract]
import { execFileSync } from 'node:child_process';
import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { dirname, join, relative } from 'node:path';
import ts from '@typescript/typescript6';

// The repository's root (git's) and Seq*'s folder in it: the scripts work wherever packages/sequences sits, its own
// repository included.
export const SEQSTAR = join(import.meta.dirname, '..');
export const root = execFileSync('git', ['rev-parse', '--show-toplevel'], { cwd: SEQSTAR, encoding: 'utf8' }).trim();
/** Seq*'s folder as git paths start: "packages/sequences/". */
export const PREFIX = relative(root, SEQSTAR) ? `${relative(root, SEQSTAR)}/` : '';

export interface Entry {
  readonly pkg: string;
  readonly subpath: string;
  readonly file: string;
}

/** An `exports` value: a path, a conditions object (the source is `molstar-src`'s) or null (hidden). */
type ExportTarget = string | Record<string, string> | null;

/** Every `.ts` entry of every package of packages/sequences, a `*` subpath expanded to the files that match it. */
export function entries(dir = SEQSTAR): Entry[] {
  const out: Entry[] = [];
  for (const pkg of readdirSync(dir)) {
    const json = join(dir, pkg, 'package.json');
    if (!existsSync(json)) continue;
    const map = JSON.parse(readFileSync(json, 'utf8')).exports as Record<string, ExportTarget> | undefined;
    if (!map) continue;
    for (const [subpath, value] of Object.entries(map)) {
      // A null target hides a subpath (`./_test/*`, `./_stories/*`, as Mol* v6); a conditions object names its source
      // under `molstar-src`. The core's `algorithms` is the vendored aligners (VENDORED.md): its surface is its own.
      const target = typeof value === 'string' ? value : value?.['molstar-src'];
      if (!target || !target.endsWith('.ts') || (pkg === 'core' && subpath.startsWith('./algorithms'))) continue;
      if (!target.includes('*')) {
        out.push({ pkg, subpath, file: join(dir, pkg, target) });
        continue;
      }
      const [head, tail] = target.split('*');
      const folder = join(dir, pkg, head);
      if (!existsSync(folder)) continue;
      for (const f of readdirSync(folder).filter((f) => f.endsWith(tail) && !f.endsWith('.d.ts'))) {
        const sub = subpath.replace('*', f.slice(0, -tail.length));
        // A file the map hides by name (`"./renderers/rows-tiles": null`) is not an entry.
        if (map[sub] !== null) out.push({ pkg, subpath: sub, file: join(folder, f) });
      }
    }
  }
  return out;
}

export interface Exported {
  readonly name: string;
  readonly space: 'type' | 'value';
  /** Where the name is declared: every declaration, `file:position`, sorted. */
  readonly decl: string;
  readonly entry: Entry;
}

const OPTIONS: ts.CompilerOptions = {
  target: ts.ScriptTarget.ES2022,
  module: ts.ModuleKind.ESNext,
  moduleResolution: ts.ModuleResolutionKind.Bundler,
  customConditions: ['molstar-src'],
  allowImportingTsExtensions: true,
  jsx: ts.JsxEmit.ReactJSX,
  noEmit: true,
  skipLibCheck: true,
  types: [],
};

/** The exports of each entry file, per namespace, resolved to their declarations. */
export function exported(list: readonly Entry[]): Exported[] {
  const program = ts.createProgram(
    list.map((e) => e.file),
    OPTIONS,
  );
  const checker = program.getTypeChecker();
  const out: Exported[] = [];
  for (const entry of list) {
    const file = program.getSourceFile(entry.file);
    const module = file && checker.getSymbolAtLocation(file);
    if (!module) continue;
    for (const s of checker.getExportsOfModule(module)) {
      const r = s.flags & ts.SymbolFlags.Alias ? checker.getAliasedSymbol(s) : s;
      const decl = (r.declarations ?? [])
        .map((d) => `${relative(root, d.getSourceFile().fileName)}:${d.pos}`)
        .sort()
        .join(',');
      if (r.flags & ts.SymbolFlags.Type) out.push({ name: s.name, space: 'type', decl, entry });
      if (r.flags & ts.SymbolFlags.Value) out.push({ name: s.name, space: 'value', decl, entry });
    }
  }
  return out;
}

/**
 * The two UIs of a package (`./react` and `./dom`) implement one slot contract under one set of names (`ViewTabs`,
 * `PartProps`, `IssueBadge`, …; ui.md §3, §7): a host imports one of them, so their overlap is the design, not a clash.
 * Either of them against any other entry is checked as for all.
 */
const alternatives = (a: Entry, b: Entry) =>
  a.pkg === b.pkg && ['./react', './dom'].includes(a.subpath) && ['./react', './dom'].includes(b.subpath);

/** One message per name that two entries export as different declarations in one namespace (the two UIs aside). */
export function collisions(list: readonly Exported[]): string[] {
  const seen = new Map<string, Exported[]>();
  for (const e of list) seen.set(`${e.space} ${e.name}`, [...(seen.get(`${e.space} ${e.name}`) ?? []), e]);
  const out: string[] = [];
  for (const [key, all] of seen) {
    const clash = all.filter((a) => all.some((b) => a.decl !== b.decl && !alternatives(a.entry, b.entry)));
    if (clash.length)
      out.push(`${key}: ${clash.map((e) => `@molstar/seq-${e.entry.pkg}${e.entry.subpath.slice(1)}`).join(', ')}`);
  }
  return out;
}

/**
 * The contract's gaps: per export of each entry (`/testing` aside), the Seq* types its declaration's signature names
 * that no entry exports, as `entry: name mentions Type (file)`.
 */
export function contractGaps(list: readonly Entry[]): string[] {
  const program = ts.createProgram(
    list.map((e) => e.file),
    OPTIONS,
  );
  const checker = program.getTypeChecker();
  const resolve = (s: ts.Symbol) => (s.flags & ts.SymbolFlags.Alias ? checker.getAliasedSymbol(s) : s);
  const shown = new Set<ts.Symbol>();
  const by: [string, ts.Symbol[]][] = [];
  for (const e of list) {
    const file = program.getSourceFile(e.file);
    const module = file && checker.getSymbolAtLocation(file);
    if (!module) continue;
    const names = checker.getExportsOfModule(module);
    for (const s of names) shown.add(resolve(s));
    if (!e.subpath.startsWith('./testing')) by.push([`@molstar/seq-${e.pkg}${e.subpath.slice(1)}`, names]);
  }
  const hidden = (n: ts.Node) =>
    ((ts.isPropertyDeclaration(n) || ts.isMethodDeclaration(n) || ts.isGetAccessor(n)) &&
      (ts.isPrivateIdentifier(n.name) ||
        n.modifiers?.some(
          (m) => m.kind === ts.SyntaxKind.PrivateKeyword || m.kind === ts.SyntaxKind.ProtectedKeyword,
        ))) ||
    ts.isBlock(n) ||
    (ts.isPropertyDeclaration(n) && !!n.initializer && !n.type);
  const out = new Set<string>();
  for (const [entry, names] of by)
    for (const s of names)
      for (const d of resolve(s).declarations ?? []) {
        const visit = (n: ts.Node) => {
          if (hidden(n)) return;
          if (ts.isTypeReferenceNode(n) || ts.isExpressionWithTypeArguments(n)) {
            const id = ts.isTypeReferenceNode(n) ? n.typeName : n.expression;
            const left = ts.isQualifiedName(id) ? checker.getSymbolAtLocation(id.left) : undefined;
            const sym = checker.getSymbolAtLocation(ts.isQualifiedName(id) ? id.right : id);
            const t = sym && resolve(sym);
            const file = t?.declarations?.[0]?.getSourceFile().fileName ?? '';
            if (
              t &&
              !(left && shown.has(resolve(left))) &&
              file.startsWith(SEQSTAR) &&
              !file.includes('/_test/') &&
              !shown.has(t) &&
              !(t.flags & ts.SymbolFlags.TypeParameter)
            )
              out.add(`${entry}: ${s.name} mentions ${t.name} (${relative(SEQSTAR, file)})`);
          }
          ts.forEachChild(n, visit);
        };
        // A variable: its type annotation only; a function: its signature, not its body.
        if (ts.isVariableDeclaration(d)) {
          if (d.type) visit(d.type);
        } else if (ts.isFunctionDeclaration(d)) {
          for (const p of d.parameters) visit(p);
          if (d.type) visit(d.type);
          for (const p of d.typeParameters ?? []) visit(p);
        } else visit(d);
      }
  return [...out];
}

/** Per entry, the exported names (either namespace) that no consumer file mentions. */
export function unused(list: readonly Exported[]): Map<string, string[]> {
  const files = execFileSync('git', ['ls-files', '-co', '--exclude-standard', '--', '*.ts', '*.tsx'], { cwd: root })
    .toString()
    .split('\n')
    .filter((f) => f && !(f.startsWith(PREFIX) && /^[^/]+\/_test\//.test(f.slice(PREFIX.length))));
  const words = new Map<string, Set<string>>();
  for (const f of files) {
    if (!existsSync(join(root, f))) continue;
    words.set(f, new Set(readFileSync(join(root, f), 'utf8').match(/[A-Za-z_$][\w$]*/g)));
  }
  const out = new Map<string, string[]>();
  for (const e of list) {
    // A consumer is outside the entry's own folder (the core's folders were packages of their own).
    const own = `${relative(root, dirname(e.entry.file))}/`;
    const used = [...words].some(([f, w]) => !f.startsWith(own) && w.has(e.name));
    const key = `@molstar/seq-${e.entry.pkg}${e.entry.subpath.slice(1)}`;
    if (!used && !out.get(key)?.includes(e.name)) out.set(key, [...(out.get(key) ?? []), e.name]);
  }
  return out;
}

if (import.meta.main) {
  const list = exported(entries());
  if (process.argv.includes('--unused')) {
    for (const [entry, names] of unused(list)) console.log(`${entry}: ${names.join(', ')}`);
  } else if (process.argv.includes('--contract')) {
    for (const g of contractGaps(entries())) console.log(g);
  } else {
    const found = collisions(list);
    for (const f of found) console.error(`export collision, ${f}`);
    const gaps = contractGaps(entries());
    for (const g of gaps) console.error(`a public signature names an unexported type, ${g}`);
    process.exit(found.length || gaps.length ? 1 : 0);
  }
}
