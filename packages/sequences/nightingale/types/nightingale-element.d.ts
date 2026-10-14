/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the type checker sees for every `@nightingale-elements/*` package (tsconfig `paths`): each element's class, all
// mount.ts uses of them (it defines the custom elements). The packages' own declarations use extensionless relative
// imports, which NodeNext rejects under the repository's full declaration check (skipLibCheck: false); bundlers and
// Vitest still resolve the real packages (they ignore `paths`).
declare const element: CustomElementConstructor;
export default element;
