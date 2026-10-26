/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { getGlslProgram, getGlslProgramNames, registerGlslProgram } from '../webgl/programs.js';

const Name = 'glsl-programs-test';

describe('GLSL programs', () => {
  it('normalizes the sources on first use, once, and again after a new registration', () => {
    registerGlslProgram(Name, '  // a comment\nvoid main() {}\n\n\n', '/* block */\nvoid main() {}\n');
    const program = getGlslProgram(Name);
    // as before the sources were normalized lazily: a comment line becomes a line break, blank lines one
    expect(program.vert).toBe('\nvoid main() {}');
    expect(program.frag).toBe('\nvoid main() {}');
    expect(getGlslProgram(Name)).toBe(program);

    registerGlslProgram(Name, 'void main() { gl_Position = vec4(0.0); }', 'void main() {}');
    expect(getGlslProgram(Name).vert).toBe('void main() { gl_Position = vec4(0.0); }');
  });

  it('rejects `#include` when the program is used', () => {
    registerGlslProgram(Name, '#include common\nvoid main() {}', 'void main() {}');
    expect(() => getGlslProgram(Name)).toThrow(/'#include' is not supported/);
  });

  it('normalizes every built-in program as before', () => {
    for (const name of getGlslProgramNames()) {
      if (name === Name) continue;
      const { vert, frag } = getGlslProgram(name);
      for (const source of [vert, frag]) {
        expect(source).not.toMatch(/\/\//);
        expect(source).not.toMatch(/\n{2,}/);
      }
    }
  });
});
