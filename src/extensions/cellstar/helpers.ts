import { ParamDefinition } from '../../mol-util/param-definition';


/** Split entry ID (e.g. 'emd-1832') into source ('emdb') and number ('1832') */
export function splitEntryId(entryId: string) {
    const PREFIX_TO_SOURCE: { [prefix: string]: string } = { 'empiar': 'empiar', 'emd': 'emdb' };
    const [prefix, entry] = entryId.split('-');
    return {
        source: PREFIX_TO_SOURCE[prefix],
        entryNumber: entry
    };
}

/** Create entry ID (e.g. 'emd-1832') for a combination of source ('emdb') and number ('1832') */
export function createEntryId(source: string, entryNumber: string | number) {
    const SOURCE_TO_PREFIX: { [prefix: string]: string } = { 'empiar': 'empiar', 'emdb': 'emd' };
    return `${SOURCE_TO_PREFIX[source]}-${entryNumber}`;
}



/**
 * Represents a set of values to choose from, with a default value. Example:
 * ```
 * export const MyChoice = new Choice({ yes: 'I agree', no: 'Nope' }, 'yes');
 * export type MyChoiceType = Choice.Values<typeof MyChoice>; // 'yes'|'no'
 * ```
 */
 export class Choice<T extends string, D extends T> {
    readonly defaultValue: D;
    readonly options: [T, string][];
    private readonly nameDict: { [value in T]: string };
    constructor(opts: { [value in T]: string }, defaultValue: D) {
        this.defaultValue = defaultValue;
        this.options = Object.keys(opts).map(k => [k as T, opts[k as T]]);
        this.nameDict = opts;
    }
    PDSelect(defaultValue?: T, info?: ParamDefinition.Info): ParamDefinition.Select<T> {
        return ParamDefinition.Select<T>(defaultValue ?? this.defaultValue, this.options, info);
    }
    prettyName(value: T): string {
        return this.nameDict[value];
    }
}
export namespace Choice {
    export type Values<T extends Choice<any, any>> = T extends Choice<infer R, any> ? R : any;
}
