/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Strict, language-agnostic shader preprocessor for the subset the shader sources use:
 * - `#if`/`#elif` with `defined(X)`/`defined X`, `!`, `&&`, `||`, parentheses and integer
 *   comparisons (`==`, `!=`, `<`, `>`, `<=`, `>=`) between integer literals and defines with
 *   integer values; `#ifdef`, `#ifndef`, `#else`, `#endif`
 * - object- and function-like `#define` (no `#`/`##` operators), `#undef`
 * - `#pragma unroll_loop_start`/`unroll_loop_end` around `for (int i = N; i < M; ++i) { ... }`
 *
 * Anything else throws: unknown directives or pragmas, other operators, undefined names or
 * non-integer values in expressions, unbalanced conditionals.
 */

/** Initial defines of a program: name to replacement text (`''` for flags) */
export type PreprocessorDefines = ReadonlyMap<string, string>;

export interface PreprocessOptions {
  /**
   * Expand macros in the code and drop `#define`/`#undef` (no preprocessor downstream, e.g. WGSL);
   * otherwise the active `#define`/`#undef` lines are kept for the GLSL compiler
   */
  readonly expandMacros?: boolean;
}

export class ShaderPreprocessorError extends Error {
  constructor(message: string, line?: number) {
    super(line === undefined ? message : `${message} (line ${line})`);
    this.name = 'ShaderPreprocessorError';
  }
}

interface Macro {
  readonly params?: readonly string[];
  readonly body: string;
}

type Macros = Map<string, Macro>;

interface Conditional {
  /** this branch is active (all enclosing conditionals are active too) */
  active: boolean;
  /** a branch of this conditional was taken */
  taken: boolean;
  /** enclosing conditionals are active */
  readonly parentActive: boolean;
  seenElse: boolean;
  readonly line: number;
}

const reDirective = /^\s*#\s*(\w+)\s*(.*)$/;
const reDefine = /^([A-Za-z_]\w*)(\(([^)]*)\))?(?:\s+(.*))?$/;
const reIdentifier = /^[A-Za-z_]\w*$/;
const reInteger = /^\d+$/;

/**
 * Resolves conditionals (and optionally macros) of `source` for the given `defines`.
 * Line breaks inside directives (`\`) are not supported.
 */
export function preprocessShader(
  source: string,
  defines: PreprocessorDefines,
  options: PreprocessOptions = {},
): string {
  const expandMacros = !!options.expandMacros;
  const macros: Macros = new Map();
  for (const [name, value] of defines) macros.set(name, { body: value });

  const stack: Conditional[] = [];
  const isActive = () => stack.length === 0 || stack[stack.length - 1].active;

  const out: string[] = [];
  const lines = unrollLoops(source).split('\n');
  for (let i = 0; i < lines.length; ++i) {
    const line = lines[i];
    const lineNumber = i + 1;
    const m = line.match(reDirective);
    if (!m) {
      if (isActive()) out.push(expandMacros ? expand(line, macros, new Set(), lineNumber) : line);
      continue;
    }

    const directive = m[1];
    const rest = m[2].replace(/\/\/.*$/, '').trim();
    const top = stack[stack.length - 1];
    switch (directive) {
      case 'if':
      case 'ifdef':
      case 'ifndef': {
        const parentActive = isActive();
        let value = false;
        if (parentActive) {
          if (directive === 'if') {
            value = evaluate(rest, macros, lineNumber);
          } else {
            if (!reIdentifier.test(rest)) throw new ShaderPreprocessorError(`#${directive} expects a name`, lineNumber);
            value = directive === 'ifdef' ? macros.has(rest) : !macros.has(rest);
          }
        }
        stack.push({ active: value, taken: value, parentActive, seenElse: false, line: lineNumber });
        break;
      }
      case 'elif': {
        if (!top) throw new ShaderPreprocessorError('#elif without #if', lineNumber);
        if (top.seenElse) throw new ShaderPreprocessorError('#elif after #else', lineNumber);
        const value = top.parentActive && !top.taken && evaluate(rest, macros, lineNumber);
        top.active = value;
        top.taken ||= value;
        break;
      }
      case 'else': {
        if (!top) throw new ShaderPreprocessorError('#else without #if', lineNumber);
        if (top.seenElse) throw new ShaderPreprocessorError('#else after #else', lineNumber);
        if (rest) throw new ShaderPreprocessorError('unexpected text after #else', lineNumber);
        top.active = top.parentActive && !top.taken;
        top.taken = true;
        top.seenElse = true;
        break;
      }
      case 'endif': {
        if (!top) throw new ShaderPreprocessorError('#endif without #if', lineNumber);
        if (rest) throw new ShaderPreprocessorError('unexpected text after #endif', lineNumber);
        stack.pop();
        break;
      }
      case 'define': {
        if (!isActive()) break;
        const d = rest.match(reDefine);
        if (!d) throw new ShaderPreprocessorError(`invalid #define '${rest}'`, lineNumber);
        const [, name, , params, body = ''] = d;
        if (body.includes('#'))
          throw new ShaderPreprocessorError(`'#' and '##' are not supported in macros`, lineNumber);
        const macro: Macro =
          params === undefined
            ? { body: body.trim() }
            : {
                params: params
                  .split(',')
                  .map((p) => p.trim())
                  .filter((p) => p.length > 0),
                body: body.trim(),
              };
        if (macro.params?.some((p) => !reIdentifier.test(p))) {
          throw new ShaderPreprocessorError(`invalid parameters of macro '${name}'`, lineNumber);
        }
        macros.set(name, macro);
        if (!expandMacros) out.push(line);
        break;
      }
      case 'undef': {
        if (!isActive()) break;
        if (!reIdentifier.test(rest)) throw new ShaderPreprocessorError('#undef expects a name', lineNumber);
        macros.delete(rest);
        if (!expandMacros) out.push(line);
        break;
      }
      case 'pragma':
        // `unroll_loop_*` pragmas were consumed by `unrollLoops`
        throw new ShaderPreprocessorError(`unsupported pragma '${rest}'`, lineNumber);
      default:
        throw new ShaderPreprocessorError(`unsupported directive '#${directive}'`, lineNumber);
    }
  }
  if (stack.length > 0) {
    throw new ShaderPreprocessorError('missing #endif', stack[stack.length - 1].line);
  }
  return out.join('\n');
}

//

const reUnrollLoop =
  /#pragma unroll_loop_start\s+for\s*\(\s*int\s+i\s*=\s*(\d+)\s*;\s*i\s*<\s*(\d+)\s*;\s*\+\+i\s*\s*\)\s*{([\s\S]+?)}\s+#pragma unroll_loop_end/g;

/**
 * Unrolls `for (int i = N; i < M; ++i) { ... }` loops (integer literal bounds) between
 * `#pragma unroll_loop_start` and `#pragma unroll_loop_end`; `[i]` and `UNROLLED_LOOP_INDEX`
 * in the body become the iteration index
 */
export function unrollLoops(str: string) {
  return str.replace(reUnrollLoop, loopReplacer);
}

function loopReplacer(_match: string, start: string, end: string, snippet: string) {
  let out = '';
  for (let i = parseInt(start); i < parseInt(end); ++i) {
    out += snippet.replace(/\[\s*i\s*\]/g, `[${i}]`).replace(/UNROLLED_LOOP_INDEX/g, `${i}`);
  }
  return out;
}

//

const reToken = /\s*(defined\b|[A-Za-z_]\w*|\d+(?:\.\d*)?|&&|\|\||==|!=|<=|>=|[<>!()])/y;

const BinaryPrecedence: { [op: string]: number } = {
  '||': 1,
  '&&': 2,
  '==': 3,
  '!=': 3,
  '<': 4,
  '>': 4,
  '<=': 4,
  '>=': 4,
};

function tokenize(expression: string, line: number) {
  const tokens: string[] = [];
  reToken.lastIndex = 0;
  let index = 0;
  while (index < expression.length) {
    if (/^\s*$/.test(expression.slice(index))) break;
    reToken.lastIndex = index;
    const m = reToken.exec(expression);
    if (!m)
      throw new ShaderPreprocessorError(
        `unsupported token in '${expression}' at '${expression.slice(index).trim()}'`,
        line,
      );
    tokens.push(m[1]);
    index = reToken.lastIndex;
  }
  return tokens;
}

/** Evaluates a `#if`/`#elif` expression of the supported subset */
function evaluate(expression: string, macros: Macros, line: number): boolean {
  if (!expression) throw new ShaderPreprocessorError('empty expression', line);
  const tokens = tokenize(expression, line);
  let i = 0;
  const fail = (message: string): never => {
    throw new ShaderPreprocessorError(`${message} in '${expression}'`, line);
  };
  const next = () => (i < tokens.length ? tokens[i++] : fail('unexpected end'));

  const integerValue = (name: string): number => {
    const macro = macros.get(name);
    if (!macro) return fail(`undefined name '${name}'`);
    if (macro.params) return fail(`function-like macro '${name}'`);
    if (!reInteger.test(macro.body)) return fail(`'${name}' is not an integer ('${macro.body}')`);
    return parseInt(macro.body);
  };

  const unary = (): number => {
    const t = next();
    if (t === '!') return unary() ? 0 : 1;
    if (t === '(') {
      const v = binary(0);
      if (next() !== ')') fail(`expected ')'`);
      return v;
    }
    if (t === 'defined') {
      const parens = tokens[i] === '(';
      if (parens) i++;
      const name = next();
      if (!reIdentifier.test(name)) fail('defined expects a name');
      if (parens && next() !== ')') fail(`expected ')'`);
      return macros.has(name) ? 1 : 0;
    }
    if (reInteger.test(t)) return parseInt(t);
    if (reIdentifier.test(t)) return integerValue(t);
    return fail(`unexpected '${t}'`);
  };

  const binary = (minPrecedence: number): number => {
    let left = unary();
    while (i < tokens.length && (BinaryPrecedence[tokens[i]] ?? 0) > minPrecedence) {
      const op = tokens[i++];
      const right = binary(BinaryPrecedence[op]);
      left = applyBinary(op, left, right);
    }
    return left;
  };

  const value = binary(0);
  if (i < tokens.length) fail(`unexpected '${tokens[i]}'`);
  return value !== 0;
}

function applyBinary(op: string, a: number, b: number): number {
  switch (op) {
    case '||':
      return a || b ? 1 : 0;
    case '&&':
      return a && b ? 1 : 0;
    case '==':
      return a === b ? 1 : 0;
    case '!=':
      return a !== b ? 1 : 0;
    case '<':
      return a < b ? 1 : 0;
    case '>':
      return a > b ? 1 : 0;
    case '<=':
      return a <= b ? 1 : 0;
    case '>=':
      return a >= b ? 1 : 0;
  }
  throw new ShaderPreprocessorError(`unsupported operator '${op}'`);
}

//

const reWord = /[A-Za-z_]\w*/g;

/** Expands macros in `text`; `disabled` holds the macros being expanded (no recursion) */
function expand(text: string, macros: Macros, disabled: ReadonlySet<string>, line: number): string {
  let out = '';
  let index = 0;
  reWord.lastIndex = 0;
  for (;;) {
    reWord.lastIndex = index;
    const m = reWord.exec(text);
    if (!m) break;
    const name = m[0];
    const start = m.index;
    let end = start + name.length;
    // skip numbers with suffixes like `1e5` or `0x1F` (the word starts inside a number)
    if (start > 0 && /[\w.]/.test(text[start - 1])) {
      out += text.slice(index, end);
      index = end;
      continue;
    }
    const macro = macros.get(name);
    if (!macro || disabled.has(name)) {
      out += text.slice(index, end);
      index = end;
      continue;
    }
    const inner = new Set(disabled).add(name);
    if (!macro.params) {
      out += text.slice(index, start) + expand(macro.body, macros, inner, line);
      index = end;
      continue;
    }
    // function-like: only an invocation if followed by `(`
    const open = text.slice(end).match(/^\s*\(/);
    if (!open) {
      out += text.slice(index, end);
      index = end;
      continue;
    }
    const { args, close } = parseArguments(text, end + open[0].length, name, line);
    if (args.length !== macro.params.length && !(macro.params.length === 0 && args.length === 1 && !args[0])) {
      throw new ShaderPreprocessorError(
        `macro '${name}' expects ${macro.params.length} arguments, got ${args.length}`,
        line,
      );
    }
    // arguments are fully expanded before substitution, like the C preprocessor
    const expandedArgs = args.map((a) => expand(a, macros, disabled, line).trim());
    const params = macro.params;
    const substituted = macro.body.replace(reWord, (word) => {
      const p = params.indexOf(word);
      return p >= 0 ? expandedArgs[p] : word;
    });
    out += text.slice(index, start) + expand(substituted, macros, inner, line);
    index = close + 1;
    end = index;
  }
  return out + text.slice(index);
}

/** Splits the arguments of a macro invocation starting after its `(`, returns the index of the closing `)` */
function parseArguments(text: string, from: number, name: string, line: number) {
  const args: string[] = [];
  let depth = 0;
  let current = '';
  for (let i = from; i < text.length; ++i) {
    const c = text[i];
    if (c === '(') depth++;
    else if (c === ')') {
      if (depth === 0) {
        args.push(current);
        return { args, close: i };
      }
      depth--;
    } else if (c === ',' && depth === 0) {
      args.push(current);
      current = '';
      continue;
    }
    current += c;
  }
  throw new ShaderPreprocessorError(`unterminated invocation of macro '${name}'`, line);
}
