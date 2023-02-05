/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 *
 * Manage dependencies which are not listed in `package.json` for performance reasons.
 */


const _loadedExtraPackages: { [dep: string]: any } = {};

/** Define imports that only get imported when first needed.
 * Example usage:
 * ```
 * const lazyImports = LazyImports.create('gl', 'jpeg-js', 'pngjs') as {
 *     'gl': typeof import('gl'),
 *     'jpeg-js': typeof import('jpeg-js'),
 *     'pngjs': typeof import('pngjs'),
 * };
 * ...
 * lazyImports.pngjs.blablabla("I'm being imported now");
 * lazyImports.pngjs.blablabla("I'm cached :D");
 * ```
 */
export class LazyImports {
    static create<U extends string>(...packages: U[]): { [dep in U]: any } {
        return new LazyImports(packages) as any;
    }
    private constructor(private packages: string[]) {
        for (const p of packages) {
            Object.defineProperty(this, p, {
                get: () => this.getPackage(p),
            });
        }
    }
    private getPackage(packageName: string) {
        if (!_loadedExtraPackages[packageName]) {
            try {
                _loadedExtraPackages[packageName] = require(packageName);
            } catch {
                const message = `Package '${packageName}' is not installed. (Some packages are not listed in the 'molstar' package dependencies for performance reasons. If you're seeing this error, you'll probably need them. If your project depends on 'molstar', add these to your dependencies: ${this.packages.join(', ')}. If you're running 'molstar' directly, run this: npm install --no-save ${this.packages.join(' ')})`;
                throw new Error(message);
            }
        }
        return _loadedExtraPackages[packageName];
    }
}
