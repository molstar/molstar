import { ParamDefinition as PD } from './molstar-lib-imports';


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
    PDSelect(defaultValue?: T, info?: PD.Info): PD.Select<T> {
        return PD.Select<T>(defaultValue ?? this.defaultValue, this.options, info);
    }
    prettyName(value: T): string {
        return this.nameDict[value];
    }
}
export namespace Choice {
    export type Values<T extends Choice<any, any>> = T extends Choice<infer R, any> ? R : any;
}
