/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Seq* spec check (packages/sequences/_spec/README.md, "Citing the spec"), run before every commit:
//   - every `seqspec:<file>#<slug>` in packages/sequences (code, tests, Storybook, the spec itself) names a heading of
//     packages/sequences/_spec/<file>.md, the slug being GitHub's anchor of the heading;
//   - the backticked identifiers and issue codes of the spec's prose, and the names declared in its `ts` code blocks,
//     are found in some Seq* package's src (the generated references in `_spec/reference/` are the tests' job);
//   - every `.ts`, `.tsx` and `.css` source file starts with Mol*'s header with the Anyfold line (AGENTS.md, Headers),
//     but the vendored, copied and generated ones (HEADER_SKIP), which keep their origin's;
//   - every issue code a package's src names (`svs.*`, `seqv.*`, `seqql.*`) is in a generated reference's code table,
//     and one written with its severity (`Issue.warning("svs.x", …)`, `{ code: "svs.x", severity: "warning" }`) has
//     the registered one (seqspec:format#validation: a code has one severity, the registry's);
//   - every heap guard file (`*.heap.ts`) runs the guard template's positive control (`controlCaught`, called in its
//     code, AGENTS.md, Tests), but those `guards-pending.txt` lists, not converted yet; a listed file that calls it is
//     a miss too;
//   - every Node heap guard file (not under `canvas/`, `dom/` or `react/`: the browser's stay deep-only) has its quick
//     twin beside it (`x.heap.ts`, `x.quick.ts`), and every `*.quick.ts` calls `quickCase` (core/_test/model/quick.ts)
//     in its code, each `control: false` with its `why` (AGENTS.md, Tests).
// Exit 1 on a miss, printing each with where it is.
//
// Usage: node packages/sequences/scripts/spec.ts check
import { execFileSync } from 'node:child_process';
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs';
import { join, relative } from 'node:path';

const SEQSTAR = join(import.meta.dirname, '..');
const root = execFileSync('git', ['rev-parse', '--show-toplevel'], { cwd: SEQSTAR, encoding: 'utf8' }).trim();
const SPEC = join(SEQSTAR, '_spec');
const CITE = /seqspec:([A-Za-z-]+)#([a-z0-9_-]+)/g;
const NAMESPACES = /^(svs|model|seqv|interactions|molstar|nightingale)\.[a-z_]+(\.[a-z_]+)+$/;
const IDENT = /^[A-Za-z_$][\w$]*(\.[A-Za-z_$][\w$]*)*$/;
const STOP = new Set(['null', 'true', 'false', 'undefined', 'number', 'string', 'boolean', 'object', 'unknown']);
// Names the spec uses that are not Seq*'s own: file formats' fields and other libraries' API.
const EXTERNAL = new Set(['entity_poly_seq']);

function walk(dir: string, keep: (f: string) => boolean, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    if (name === 'node_modules' || name === 'dist' || name === 'lib') continue;
    const f = join(dir, name);
    if (statSync(f).isDirectory()) walk(f, keep, out);
    else if (keep(f)) out.push(f);
  }
  return out;
}

/** GitHub's anchor of a heading: lowercased, markup and punctuation dropped, spaces to hyphens. */
export function slugOf(heading: string): string {
  return heading
    .replace(/\\(.)/g, '$1')
    .replace(/`/g, '')
    .trim()
    .toLowerCase()
    .replace(/[^\p{L}\p{N}\s_-]/gu, '')
    .replace(/\s/g, '-');
}

/** Each spec file's anchors, a repeated heading's as GitHub numbers them (`-1`, `-2`). */
function anchors(): Map<string, Set<string>> {
  const out = new Map<string, Set<string>>();
  for (const name of readdirSync(SPEC).filter((f) => f.endsWith('.md'))) {
    const set = new Set<string>();
    let fence = false;
    for (const line of readFileSync(join(SPEC, name), 'utf8').split('\n')) {
      if (line.startsWith('```')) fence = !fence;
      const m = !fence && /^#{1,6}\s+(.*)$/.exec(line);
      if (!m) continue;
      const base = slugOf(m[1]);
      let slug = base;
      for (let n = 1; set.has(slug); n++) slug = `${base}-${n}`;
      set.add(slug);
    }
    out.set(name.slice(0, -3), set);
  }
  return out;
}

function citations(misses: string[]): number {
  const heads = anchors();
  let count = 0;
  for (const file of walk(SEQSTAR, (f) => /\.(ts|tsx|mdx|md|css|grit|json)$/.test(f))) {
    const lines = readFileSync(file, 'utf8').split('\n');
    for (let i = 0; i < lines.length; i++) {
      for (const [cite, stem, slug] of lines[i].matchAll(CITE)) {
        count++;
        const set = heads.get(stem);
        if (!set) misses.push(`${relative(root, file)}:${i + 1}\t${cite}: no spec file _spec/${stem}.md`);
        else if (!set.has(slug)) misses.push(`${relative(root, file)}:${i + 1}\t${cite}: no such heading`);
      }
    }
  }
  return count;
}

/** The names a line of TypeScript declares: types, functions, constants and members. */
function declared(text: string): string[] {
  const code = text.replace(/\/\/.*$/, '').replace(/\/\*.*?\*\//g, '');
  const out: string[] = [];
  for (const m of code.matchAll(/\b(?:interface|type|const|function|namespace|class)\s+([A-Za-z_$][\w$]*)/g))
    out.push(m[1]);
  for (const m of code.matchAll(/(?:^|[{;,(]|\breadonly)\s*([A-Za-z_$][\w$]*)\??\s*(?:<[^>]*>)?\s*[:(]/g))
    out.push(m[1]);
  return out.filter((n) => !STOP.has(n) && n.length > 2 && !/^(readonly|declare)$/.test(n));
}

/** What to look for in one backticked span: an issue code, or the names of an identifier; nothing for values. */
function namesOf(span: string): string[] {
  const s = span
    .trim()
    .replace(/\(.*\)$/s, '')
    .replace(/\[.*\]$/s, '')
    .replace(/:.*$/s, '')
    .trim();
  if (NAMESPACES.test(s)) return [s];
  if (!IDENT.test(s) || STOP.has(s) || EXTERNAL.has(s)) return [];
  if (/^[A-Z]+$/.test(s) || s.length <= 2) return []; // residue strings, single letters
  return [s];
}

/** The Mol* adapter (libraries.md, The Mol* adapter) is Prism's: its names count as Seq*'s. */
const ADAPTER = join(SEQSTAR, '../prism/viewer/src/sequence/molstar');

function identifiers(misses: string[]): number {
  const texts = [...readdirSync(SEQSTAR).map((pkg) => join(SEQSTAR, pkg, 'src')), ADAPTER]
    .filter((src) => existsSync(src))
    .flatMap((src) => walk(src, (f) => /\.(ts|tsx|js|json|css)$/.test(f)))
    .map((f) => readFileSync(f, 'utf8'));
  const corpus = texts.join('\n');
  const word = (w: string) => new RegExp(`(^|[^\\w$])${w.replace(/\$/g, '\\$')}($|[^\\w$])`);
  const found = new Map<string, boolean>();
  const has = (name: string): boolean => {
    let r = found.get(name);
    if (r !== undefined) return r;
    if (NAMESPACES.test(name)) r = corpus.includes(name);
    else {
      // `A.b.c`: every part a word, and the last one in a file that also names the first.
      const parts = name.split('.');
      const last = word(parts[parts.length - 1]);
      const first = word(parts[0]);
      r = parts.every((p) => word(p).test(corpus)) && texts.some((t) => last.test(t) && first.test(t));
    }
    found.set(name, r);
    return r;
  };
  let count = 0;
  for (const name of readdirSync(SPEC).filter((f) => f.endsWith('.md'))) {
    const file = join(SPEC, name);
    const where = (i: number) => `${relative(root, file)}:${i + 1}`;
    const note = (n: string, i: number) => {
      count++;
      if (!has(n)) misses.push(`${where(i)}\t\`${n}\` is in no Seq* src`);
    };
    let ts = false;
    let fence = false;
    const lines = readFileSync(file, 'utf8').split('\n');
    for (let i = 0; i < lines.length; i++) {
      const text = lines[i];
      if (text.startsWith('```')) {
        fence = !fence;
        ts = fence && /^```(ts|typescript)\b/.test(text);
        continue;
      }
      if (fence) {
        if (ts) for (const n of declared(text)) note(n, i);
        continue;
      }
      // A code written short after a full one on its line (`svs.source.load`, `.format`) is that code's sibling.
      let code = '';
      for (const m of text.matchAll(/`([^`]+)`/g)) {
        const span = m[1].trim();
        if (NAMESPACES.test(span)) code = span;
        else if (code && /^\.[a-z_]+$/.test(span)) {
          note(code.slice(0, code.lastIndexOf('.')) + span, i);
          continue;
        }
        for (const n of namesOf(span)) note(n, i);
      }
    }
  }
  return count;
}

/** The code tables of the generated references: `| \`code\` | severity |`. */
function registeredCodes(): Map<string, string> {
  const out = new Map<string, string>();
  const dir = join(SPEC, 'reference');
  for (const name of readdirSync(dir).filter((f) => f.endsWith('.md'))) {
    const text = readFileSync(join(dir, name), 'utf8');
    for (const [, code, severity] of text.matchAll(/^\| `([a-z]+(?:\.[a-z_]+)+)` \| (error|warning|info) \|/gm))
      out.set(code, severity);
  }
  return out;
}

/** The namespaces whose codes the references list; other namespaces (an extension's) register their own. */
const REFERENCED = /^(svs|seqv|seqql)\./;
// One part or more after the namespace (LA-19: `seqql.type`, `seqql.empty`, … went unchecked).
const CODE = /["'`]((?:svs|seqv|seqql)\.[a-z_]+(?:\.[a-z_-]+)*)["'`]/g;
const WRITTEN = [
  /Issue\.(error|warning|info)\(\s*["']([a-z]+(?:\.[a-z_-]+)+)["']/g,
  /code: ["']([a-z]+(?:\.[a-z_-]+)+)["'],\s*severity: ["'](error|warning|info)["']/g,
];

function codes(misses: string[]): number {
  const registered = registeredCodes();
  const files = readdirSync(SEQSTAR)
    .map((pkg) => join(SEQSTAR, pkg, 'src'))
    .filter((src) => existsSync(src))
    .flatMap((src) => walk(src, (f) => /\.(ts|tsx)$/.test(f)));
  let count = 0;
  for (const file of files) {
    const text = readFileSync(file, 'utf8');
    const line = (at: number) => `${relative(root, file)}:${text.slice(0, at).split('\n').length}`;
    for (const m of text.matchAll(CODE)) {
      count++;
      if (!registered.has(m[1])) misses.push(`${line(m.index)}\t"${m[1]}" is in no reference's issue codes`);
    }
    for (const [k, re] of WRITTEN.entries()) {
      for (const m of text.matchAll(re)) {
        const [code, severity] = k === 0 ? [m[2], m[1]] : [m[1], m[2]];
        const want = registered.get(code);
        if (REFERENCED.test(code) && want && want !== severity)
          misses.push(`${line(m.index)}\t"${code}" is written as ${severity}; it is registered as ${want}`);
      }
    }
  }
  return count;
}

/**
 * A TypeScript file's text with its comments and its strings' and templates' text blanked (a template's `${…}` kept,
 * as code); regular expression literals are not told apart from division, which the guard files do not need.
 */
function codeOnly(text: string): string {
  let out = '';
  const braces: number[] = []; // per open template expression, its `{` depth
  let depth = 0;
  for (let i = 0; i < text.length; ) {
    const c = text[i];
    const two = text.slice(i, i + 2);
    if (two === '//') {
      while (i < text.length && text[i] !== '\n') i++;
    } else if (two === '/*') {
      const end = text.indexOf('*/', i + 2);
      i = end < 0 ? text.length : end + 2;
      out += ' ';
    } else if (c === '"' || c === "'") {
      for (i++; i < text.length && text[i] !== c && text[i] !== '\n'; i++) if (text[i] === '\\') i++;
      i++;
      out += '""';
    } else if (c === '`' || (c === '}' && braces.length > 0 && braces[braces.length - 1] === depth)) {
      // A template's text, from its start or the end of an expression in it, to its end or its next `${`.
      if (c === '}') braces.pop();
      for (i++; i < text.length && text[i] !== '`'; i++) {
        if (text[i] === '\\') i++;
        else if (text[i] === '$' && text[i + 1] === '{') break;
      }
      if (text[i] === '$') {
        braces.push(depth);
        i += 2;
      } else i++;
      out += '``';
    } else {
      if (c === '{') depth++;
      else if (c === '}') depth--;
      out += c;
      i++;
    }
  }
  return out;
}

/** The guard files not on the template yet, one path under packages/sequences per line (`#` comments). */
const PENDING = join(SEQSTAR, 'scripts', 'guards-pending.txt');

function guards(misses: string[]): number {
  const pending = new Set(
    readFileSync(PENDING, 'utf8')
      .split('\n')
      .map((l) => l.replace(/#.*$/, '').trim())
      .filter((l) => l.length > 0),
  );
  const files = walk(SEQSTAR, (f) => f.endsWith('.heap.ts')).map((f) => relative(SEQSTAR, f));
  for (const f of files) {
    // A call in the code, not a mention in a comment or a string (LC-23: `/* controlCaught( */` passed).
    const templated = /\bcontrolCaught\(/.test(codeOnly(readFileSync(join(SEQSTAR, f), 'utf8')));
    if (!templated && !pending.has(f))
      misses.push(
        `${f}: a heap guard without the template's control (controlCaught; packages/sequences/AGENTS.md, Tests)`,
      );
    if (templated && pending.has(f))
      misses.push(`${f}: on the template now; remove it from scripts/guards-pending.txt`);
  }
  for (const f of pending)
    if (!files.includes(f)) misses.push(`scripts/guards-pending.txt: "${f}" is no heap guard file`);
  for (const f of files) {
    const twin = f.replace(/\.heap\.ts$/, '.quick.ts');
    if (!BROWSER_GUARD.test(f) && !existsSync(join(SEQSTAR, twin)))
      misses.push(
        `${f}: a Node heap guard without its quick twin ${twin} (quickCase; packages/sequences/AGENTS.md, Tests)`,
      );
  }
  const quick = walk(SEQSTAR, (f) => f.endsWith('.quick.ts')).map((f) => relative(SEQSTAR, f));
  for (const f of quick) {
    const code = codeOnly(readFileSync(join(SEQSTAR, f), 'utf8'));
    if (!/\bquickCase\(/.test(code))
      misses.push(`${f}: a quick file that never calls quickCase (core/_test/model/quick.ts; AGENTS.md, Tests)`);
    const off = code.match(/\bcontrol:\s*false\b/g)?.length ?? 0;
    const why = code.match(/\bwhy:/g)?.length ?? 0;
    if (off > why) misses.push(`${f}: a quick case without its control says why (\`control: false, why: "…"\`)`);
  }
  return files.length + quick.length;
}

/** The browser's heap guards (`test:browser-heap`), deep-only: no `v8.GCProfiler` in Chrome. */
const BROWSER_GUARD = /(^|\/)_test\/(.*\/)?(canvas|dom|react)\//;

const HEADER = `/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
`;
/** Vendored or copied (each folder's VENDORED.md) and generated files: they keep their origin's header. */
const HEADER_SKIP = [
  /^core\/src\/algorithms\/biomsa\//,
  /^storybook\/src\/_stories\/ui\/prism\/(reset|theme)\.css$/,
  /^viewer\/src\/svg\/metrics-tables\.ts$/,
];

/** Source files that do not start with the header (after a `#!` line). */
function headers(misses: string[]): number {
  const files = walk(SEQSTAR, (f) => /\.(ts|tsx|css)$/.test(f))
    .map((f) => relative(SEQSTAR, f))
    .filter((f) => !HEADER_SKIP.some((r) => r.test(f)));
  for (const f of files) {
    const text = readFileSync(join(SEQSTAR, f), 'utf8');
    const body = text.startsWith('#!') ? text.slice(text.indexOf('\n') + 1) : text;
    if (!body.startsWith(HEADER))
      misses.push(`${f}: no Mol* header with the Anyfold line (packages/sequences/AGENTS.md)`);
  }
  return files.length;
}

function main(): void {
  const cmd = process.argv[2];
  if (cmd !== 'check') {
    console.error('Usage: node packages/sequences/scripts/spec.ts check');
    process.exitCode = 2;
    return;
  }
  const misses: string[] = [];
  const cites = citations(misses);
  const names = identifiers(misses);
  const files = headers(misses);
  const named = codes(misses);
  const guarded = guards(misses);
  for (const m of misses) console.log(m);
  console.log(
    `${cites} citation(s), ${names} spec name(s), ${files} file header(s), ${named} issue code(s), ${guarded} guard file(s) checked; ${misses.length} miss(es)`,
  );
  if (misses.length) process.exitCode = 1;
}

main();
