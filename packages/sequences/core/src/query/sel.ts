/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `sel` shorthand (seqspec:seqql#the-sel-shorthand): refs with label ranges, index slices, item ids and inputs,
// unioned by commas, expanded at bind time to base symbols. Its errors point into the string; the stored string is
// normalized (spaces removed) once it reads cleanly.
import { NUMBER } from "./read.js";
import type { Apply, Expression, MacroReport } from "./types.js";

const LABEL_SOURCE = "-?[0-9]+(?:[A-Za-z]+|\\.[0-9]+)?|[A-Za-z][A-Za-z0-9_]*";
const LABEL = new RegExp(LABEL_SOURCE, "y");
/** A numbering's name, a namespaced one included (`uniprot:P04637`: to the last `:` before the labels; CM-02). */
const NUMBERING = /[A-Za-z][A-Za-z0-9_-]*(?::[A-Za-z0-9_.-]+)*(?=\s*:)/y;
const INT = /[0-9]+/y;
const REF = /[^\s[\],#$@]+/y;
const NAME = /[a-z][a-z0-9-]*/y;
const ID = /[^\s[\],#]+/y;
const SPACE = /\s*/y;
/** A label or label range as `range` takes it (seqspec:seqql#labels-and-indices): "10-20", "45A", "-5-10", a name. */
export const LABEL_RANGE = new RegExp(`^(?:${LABEL_SOURCE})(?:-(?:${LABEL_SOURCE}))?$`);
/** What `sel` takes, its message when the argument is missing or not a string. */
export const SEL_STRING = 'sel takes the shorthand as a string: (sel "1brs/A[10-20]")';
/** An id that does not end at "," or "]": ids are one word, and labels are matched, not listed (QE1a, q21). */
const ID_WORD = (id: string, end: string) =>
  `an id is one word: ${id} ends at ${/\s/.test(end) ? "the space" : `"${end}"`}; ids are the annotation's own, separated by commas; match labels with (pick (items :of A) :where (matches (label) …))`;
const USE = 'labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")';

type Problem = Required<Pick<MacroReport, "code" | "message" | "index">> & Pick<MacroReport, "open">;
class SelError extends Error {
  readonly problem: Problem;
  constructor(problem: Problem) {
    super(problem.message);
    this.problem = problem;
  }
}
const call = (head: string, args: Record<string, Expression | Expression[]>): Apply => ({ head, args });
/** A label as the binder reads it: a number when the word is one and prints back the same, else the string. */
const labelValue = (l: string) => (NUMBER.test(l) && String(Number(l)) === l ? Number(l) : l);

/** The expansion of a shorthand string, or its first problem (indices into the string). */
export function parseSel(s: string): { expr: Expression } | { problem: Problem } {
  let i = 0;
  const fail = (message: string, index = i, open?: number, code = "seqql.sel.syntax"): never => {
    throw new SelError({ code, message, index, ...(open === undefined ? {} : { open }) });
  };
  const take = (re: RegExp): string | null => {
    re.lastIndex = i;
    const m = re.exec(s);
    if (m) i = re.lastIndex;
    return m ? m[0] : null;
  };
  const unexpected = () => fail(i < s.length ? `unexpected "${s[i]}": ${USE}` : `the shorthand ends early: ${USE}`);
  const mixed = () => fail("a bracket holds labels or #ids, not both", i, undefined, "seqql.sel.mixed");
  /** `item { "," item }`, spaces allowed after the commas. */
  const list = (item: () => void) => {
    item();
    for (take(SPACE); s[i] === ","; take(SPACE)) {
      i++;
      take(SPACE);
      item();
    }
  };

  function bracket(base: Expression): Expression {
    const open = i++;
    take(SPACE);
    const ids = s[i] === "#";
    let [idStart, idEnd] = [i, i];
    const numbering = ids ? null : take(NUMBERING);
    if (numbering) i = s.indexOf(":", i) + 1; // after the name, the separator
    take(SPACE);
    const parts: Expression[] = [];
    list(() => {
      if (ids !== (s[i] === "#")) return ids && !/[-0-9A-Za-z]/.test(s[i] ?? "") ? unexpected() : mixed();
      if (ids) {
        idStart = i++;
        parts.push(take(ID) ?? fail("an item id follows #: [#PF00069]"));
        idEnd = i;
        return;
      }
      const start = i;
      if (!take(LABEL)) fail("expected a label: 10-20, 45A");
      if (s[i] === "-" && ++i && !take(LABEL)) fail('a range is two labels joined by "-": 10-20');
      return parts.push(labelValue(s.slice(start, i)));
    });
    if (s[i] !== "]" && ids && idEnd < s.length) fail(ID_WORD(s.slice(idStart, idEnd), s[idEnd] as string), idEnd);
    if (s[i] !== "]") fail('expected "]" to close "["', i, open);
    i++;
    const args: Record<string, Expression | Expression[]> = { "0": base, "1": parts };
    if (numbering) args.numbering = numbering;
    return call(ids ? "seq.filter.ids" : "seq.modifier.range", args);
  }

  function term(): Expression {
    let base: Expression;
    if (s[i] === "$") {
      i++;
      base = call("seq.generator.input", { "0": take(NAME) ?? fail("an input is $ and its name: $selection") });
    } else base = take(REF) ?? unexpected();
    if (s[i] === "[") return bracket(base);
    if (s[i] !== "@") return base;
    const at = ++i;
    const start = take(INT);
    const stop = start === null ? null : s[i] === ":" && ++i ? take(INT) : String(Number(start) + 1);
    if (start === null || stop === null) return fail("a slice is start:stop, indices from 0: @9:20");
    // A start after the stop wraps a cyclic system's origin, known when it runs (CM-10, as `slice`): equal ends only.
    if (Number(stop) === Number(start)) fail("a slice's ends are equal, which selects nothing: @9:20", at);
    return call("seq.modifier.slice", { "0": base, "1": [Number(start), Number(stop)] });
  }

  try {
    const terms: Expression[] = [];
    take(SPACE);
    list(() => terms.push(term()));
    if (i < s.length) unexpected();
    const [first, ...rest] = terms as [Expression, ...Expression[]];
    return { expr: rest.length === 0 ? first : call("seq.combinator.merge", { "0": first, "1": rest }) };
  } catch (e) {
    if (e instanceof SelError) return { problem: e.problem };
    throw e;
  }
}

/** `sel`'s expansion: its string parsed; a problem is reported into the string and the call selects nothing. */
export function expandSel(c: Apply, report: (r: MacroReport) => void): Expression {
  const s = c.args?.["0"];
  const r = typeof s === "string" ? parseSel(s) : null;
  if (r && "expr" in r) return r.expr;
  report(r ? r.problem : { code: s === undefined ? "seqql.arg.missing" : "seqql.type", message: SEL_STRING });
  return { head: "seq.generator.empty" };
}

/** The stored form of a `sel` call: its string without spaces when it reads cleanly, else as written. */
export function normalizeSel(c: Apply): Apply {
  const s = c.args?.["0"];
  if (typeof s !== "string" || !/\s/.test(s) || !("expr" in parseSel(s))) return c;
  return { head: c.head, args: { ...c.args, "0": s.replace(/\s+/g, "") } };
}
