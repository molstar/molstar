/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { CylindersSchema } from '../renderable/cylinders.js';
import { DirectVolumeSchema } from '../renderable/direct-volume.js';
import { ImageSchema } from '../renderable/image.js';
import { LinesSchema } from '../renderable/lines.js';
import { MeshSchema } from '../renderable/mesh.js';
import { PointsSchema } from '../renderable/points.js';
import type { RenderableSchema } from '../renderable/schema.js';
import { SpheresSchema } from '../renderable/spheres.js';
import { TextSchema } from '../renderable/text.js';
import { TextureMeshSchema } from '../renderable/texture-mesh.js';
import {
  CylindersShaderCode,
  type DefineType,
  DirectVolumeShaderCode,
  ImageShaderCode,
  LinesShaderCode,
  MeshShaderCode,
  PointsShaderCode,
  type ShaderCode,
  SpheresShaderCode,
  TextShaderCode,
} from '../shader-code.js';
import { preprocessShader } from '../shader-preprocessor.js';

/** test helpers shared by the shader tests, not a test file itself */

export const GraphicsShaders: [name: string, shaderCode: ShaderCode, schema: RenderableSchema, geometryType: string][] =
  [
    ['mesh', MeshShaderCode, MeshSchema, 'mesh'],
    ['texture-mesh', MeshShaderCode, TextureMeshSchema, 'textureMesh'],
    ['spheres', SpheresShaderCode, SpheresSchema, 'spheres'],
    ['cylinders', CylindersShaderCode, CylindersSchema, 'cylinders'],
    ['points', PointsShaderCode, PointsSchema, 'points'],
    ['lines', LinesShaderCode, LinesSchema, 'lines'],
    ['text', TextShaderCode, TextSchema, 'text'],
    ['image', ImageShaderCode, ImageSchema, 'image'],
    ['direct-volume', DirectVolumeShaderCode, DirectVolumeSchema, 'directVolume'],
  ];

/**
 * Flattened define names (with values for number defines) as set by `addShaderDefines`,
 * including the shader's `ignoreDefine` filter
 */
export function getProgramDefines(
  shaderCode: ShaderCode,
  defineValues: { [k: string]: DefineType | number },
  prefixDefines: readonly string[],
): Map<string, string> {
  const cells: { [k: string]: ValueCell<any> } = {};
  for (const k of Object.keys(defineValues)) cells[k] = ValueCell.create(defineValues[k]);
  const variant = `${defineValues.dRenderVariant ?? ''}`;
  const defines = new Map<string, string>();
  for (const name of prefixDefines) defines.set(name, '');
  for (const [k, v] of Object.entries(defineValues)) {
    if (shaderCode.ignoreDefine?.(k, variant, cells)) continue;
    if (typeof v === 'string') defines.set(`${k}_${v}`, '');
    else if (typeof v === 'number') defines.set(k, `${v}`);
    else if (v) defines.set(k, '');
  }
  return defines;
}

/**
 * Preprocesses a program source like `addShaderDefines`: count defines are replaced in the
 * source (also when ignored), then the shader preprocessor resolves the conditionals;
 * returns the active lines without `#define`/`#undef`
 */
export function preprocessProgram(
  source: string,
  defines: ReadonlyMap<string, string>,
  values: { [k: string]: DefineType | number },
): string {
  for (const name of ['dLightCount']) {
    if (values[name] !== undefined) source = source.replace(new RegExp(name, 'g'), `${values[name]}`);
  }
  return preprocessShader(source, defines)
    .split('\n')
    .filter((l) => !/^\s*#\s*(define|undef)\b/.test(l))
    .join('\n');
}

/**
 * Lenient reference GLSL preprocessor for the shader preprocessor tests: evaluates conditionals
 * (any arithmetic, undefined names are 0) and tracks `#define`/`#undef`;
 * returns the active lines without directives (except `#version`/`#extension`/`#pragma`)
 */
export function preprocessGlsl(source: string, initial: ReadonlyMap<string, string>): string {
  const macros = new Map(initial);
  const stack: { active: boolean; taken: boolean; parentActive: boolean }[] = [];
  const isActive = () => stack.every((f) => f.active);
  const evaluate = (expression: string) => {
    const substituted = expression
      .replace(/defined\s*\(\s*(\w+)\s*\)|defined\s+(\w+)/g, (_, a, b) => (macros.has(a ?? b) ? '1' : '0'))
      .replace(/[A-Za-z_]\w*/g, (name) => {
        const value = macros.get(name)?.trim();
        return value && /^-?\d+(\.\d+)?$/.test(value) ? value : '0';
      });
    return evaluateExpression(substituted) !== 0;
  };

  const lines: string[] = [];
  for (const line of source.split('\n')) {
    const m = line.match(/^\s*#\s*(\w+)\s*(.*)$/);
    if (!m) {
      if (isActive()) lines.push(line);
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
        const value =
          parentActive &&
          (directive === 'ifdef' ? macros.has(rest) : directive === 'ifndef' ? !macros.has(rest) : evaluate(rest));
        stack.push({ active: value, taken: value, parentActive });
        break;
      }
      case 'elif': {
        const value = top.parentActive && !top.taken && evaluate(rest);
        top.active = value;
        top.taken ||= value;
        break;
      }
      case 'else':
        top.active = top.parentActive && !top.taken;
        top.taken = true;
        break;
      case 'endif':
        stack.pop();
        break;
      case 'define':
        if (isActive()) {
          const [, name, value] = rest.match(/^(\w+)(?:\([^)]*\))?\s*(.*)$/)!;
          macros.set(name, value);
        }
        break;
      case 'undef':
        if (isActive()) macros.delete(rest);
        break;
      default:
        if (isActive()) lines.push(line);
    }
  }
  return lines.join('\n');
}

const BinaryPrecedence: { [op: string]: number } = {
  '||': 1,
  '&&': 2,
  '==': 3,
  '!=': 3,
  '<': 4,
  '>': 4,
  '<=': 4,
  '>=': 4,
  '+': 5,
  '-': 5,
  '*': 6,
  '/': 6,
  '%': 6,
};

/** numeric preprocessor expression (identifiers already substituted), precedence climbing */
function evaluateExpression(expression: string): number {
  const tokens = expression.match(/\d+(?:\.\d+)?|&&|\|\||[=!<>]=|[-+*/%()!<>]/g) ?? [];
  let i = 0;
  const unary = (): number => {
    const t = tokens[i++];
    if (t === '!') return unary() ? 0 : 1;
    if (t === '-') return -unary();
    if (t === '+') return unary();
    if (t === '(') {
      const v = binary(0);
      i += 1; // ')'
      return v;
    }
    return Number(t);
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
  return binary(0);
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
    case '+':
      return a + b;
    case '-':
      return a - b;
    case '*':
      return a * b;
    case '/':
      return a / b;
    default:
      return a % b;
  }
}
