/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export class ErrorContext {
    private errors: { [tag: string]: any[] } = Object.create(null);

    get(tag: string): ReadonlyArray<any> {
        return this.errors[tag] ?? [];
    }

    add(tag: string, error: any) {
        if (tag in this.errors && Array.isArray(this.errors[tag])) {
            this.errors[tag].push(error);
        } else {
            this.errors[tag] = [error];
        }
    }

    clear(tag: string) {
        delete this.errors[tag];
    }
}