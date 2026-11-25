/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An `include` pattern (_spec/README.md, Language, Assemblies): `^…`, a run of atoms (a character, `.`, a class `[…]`,
// `\d`, `\w`, an escaped character), each with at most one of `?`, `*`, `+`, an optional `$` at the end. No groups,
// alternation, backreferences or counted repeats: a spec is untrusted, and a pattern is matched against every name of
// an assembly (100,000 of them), so the matcher is a table over the name and the atoms, never a backtracking search.

/** One atom of a pattern: the characters it takes, and how many times (`1`, `?`, `*`, `+`). */
interface Atom {
  readonly test: (c: string) => boolean;
  /** The one character it takes, when a plain or escaped character (a literal prefix is tested by `startsWith`). */
  readonly char?: string;
  readonly min: 0 | 1;
  readonly many: boolean;
}

/** The longest pattern taken. */
export const MaxIncludePattern = 256;

/** A class's body (`0-9XYM`, `^…` negated) as a test, or undefined when it does not read. */
function classTest(body: string): ((c: string) => boolean) | undefined {
  const neg = body.startsWith('^');
  const s = neg ? body.slice(1) : body;
  if (!s.length) return undefined;
  const ranges: [string, string][] = [];
  for (let i = 0; i < s.length; i++) {
    let a = s[i];
    if (a === '\\') a = s[++i];
    if (a === undefined) return undefined;
    if (s[i + 1] === '-' && i + 2 < s.length) {
      const b = s[i + 2];
      if (b < a) return undefined;
      ranges.push([a, b]);
      i += 2;
    } else ranges.push([a, a]);
  }
  return (c) => ranges.some(([a, b]) => c >= a && c <= b) !== neg;
}

/** The atoms of a pattern, or undefined when it is not one this language takes. */
export function includeAtoms(pattern: string): readonly Atom[] | undefined {
  if (!pattern.startsWith('^') || pattern.length > MaxIncludePattern) return undefined;
  const end = pattern.endsWith('$') && !pattern.endsWith('\\$') ? pattern.length - 1 : pattern.length;
  const atoms: Atom[] = [];
  for (let i = 1; i < end; ) {
    let test: ((c: string) => boolean) | undefined;
    let char: string | undefined;
    const c = pattern[i];
    if (c === '[') {
      const close = pattern.indexOf(']', i + 2);
      if (close < 0 || close >= end) return undefined;
      test = classTest(pattern.slice(i + 1, close));
      i = close + 1;
    } else if (c === '\\') {
      const e = pattern[i + 1];
      if (e === undefined || i + 1 >= end) return undefined;
      if (e === 'd') test = (x) => x >= '0' && x <= '9';
      else if (e === 'w') test = (x) => /\w/.test(x);
      else if (!/[a-zA-Z0-9]/.test(e)) {
        test = (x) => x === e;
        char = e;
      }
      i += 2;
    } else if ('()|{}?*+$^'.includes(c)) return undefined;
    else {
      test = c === '.' ? () => true : (x) => x === c;
      if (c !== '.') char = c;
      i++;
    }
    if (!test) return undefined;
    const q = pattern[i];
    const quantified = i < end && (q === '?' || q === '*' || q === '+');
    atoms.push({ test, char, min: quantified && q !== '+' ? 0 : 1, many: quantified && q !== '?' });
    if (quantified) i++;
  }
  return atoms;
}

/**
 * A test of names against the atoms from a name's start, to its end when `anchored`: a table of the atoms each prefix
 * of the name can reach (O(name × atoms), no backtracking), its two rows allocated once; the pattern's literal prefix
 * (`^chr`) tested first, so most names of an assembly are refused without the table.
 */
export function atomsMatcher(atoms: readonly Atom[], anchored: boolean): (name: string) => boolean {
  let lit = 0;
  while (lit < atoms.length && atoms[lit].char !== undefined && atoms[lit].min === 1 && !atoms[lit].many) lit++;
  const prefix = atoms
    .slice(0, lit)
    .map((a) => a.char)
    .join('');
  const rest = atoms.slice(lit);
  const n = rest.length;
  // `at[k]`: the name read so far is matched by the atoms before k (optional atoms skipped)
  let at = new Uint8Array(n + 1);
  let next = new Uint8Array(n + 1);
  const close = (s: Uint8Array) => {
    for (let k = 0; k < n; k++) if (s[k] && rest[k].min === 0) s[k + 1] = 1;
  };
  return (name) => {
    if (!name.startsWith(prefix)) return false;
    at.fill(0);
    at[0] = 1;
    close(at);
    if (!anchored && at[n]) return true;
    for (let i = prefix.length; i < name.length; ) {
      const c = (name.codePointAt(i) ?? 0) > 0xffff ? name.slice(i, i + 2) : name[i];
      i += c.length;
      next.fill(0);
      let any = false;
      for (let k = 0; k < n; k++)
        if (at[k] && rest[k].test(c)) {
          next[k + 1] = 1;
          if (rest[k].many) next[k] = 1;
          any = true;
        }
      // no atom took the character: no later one can
      if (!any) return false;
      close(next);
      [at, next] = [next, at];
      if (!anchored && at[n]) return true;
    }
    return at[n] === 1;
  };
}

/** Whether `name` matches the atoms (`atomsMatcher` for one name). */
export function matchAtoms(atoms: readonly Atom[], name: string, anchored: boolean): boolean {
  return atomsMatcher(atoms, anchored)(name);
}

/** An `include` pattern's test, or undefined when it is not one this language takes (`seq.lang.include.pattern`). */
export function includePattern(pattern: string): ((name: string) => boolean) | undefined {
  const atoms = includeAtoms(pattern);
  if (!atoms) return undefined;
  return atomsMatcher(atoms, pattern.endsWith('$') && !pattern.endsWith('\\$'));
}
