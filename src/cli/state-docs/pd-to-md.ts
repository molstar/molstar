/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';

export function paramsToMd(params: PD.Params) {
    return getParams(params, 0);
}

function paramInfo(param: PD.Any, offset: number): string {
    switch (param.type) {
        case 'value': return 'Value';
        case 'boolean': return 'true/false';
        case 'number': return 'Numeric value';
        case 'converted': return paramInfo(param.converted, offset);
        case 'conditioned': return getParams(param.conditionParams, offset);
        case 'multi-select': return `Array of ${oToS(param.options)}`;
        case 'color': return 'Color as 0xrrggbb';
        case 'color-list': return `A list of colors as 0xrrggbb`;
        case 'vec3': return `3D vector [x, y, z]`;
        case 'mat4': return `4x4 transformation matrix`;
        case 'url': return `URL couple with unique identifier`;
        case 'file': return `JavaScript File Handle`;
        case 'file-list': return `JavaScript FileList Handle`;
        case 'select': return `One of ${oToS(param.options)}`;
        case 'text': return 'String';
        case 'interval': return `Interval [min, max]`;
        case 'group': return `Object with:\n${getParams(param.params, offset + 2)}`;
        case 'mapped': return `Object { name: string, params: object } where name+params are:\n${getMapped(param, offset + 2)}`;
        case 'line-graph': return `A list of 2d vectors [xi, yi][]`;
        case 'object-list': return `Array of\n${paramInfo(PD.Group(param.element), offset + 2)}`;
        // TODO: support more languages
        case 'script': return `An expression in the specified language { language: 'mol-script', expressiong: string }`;
        default:
            const _: never = param;
            console.warn(`${_} has no associated UI component`);
            return '';
    }
}

function oToS(options: readonly (readonly [string, string] | readonly [string, string, string | undefined])[]) {
    return options.map(o => `'${o[0]}'`).join(', ');
}

function offsetS(n: number) {
    return new Array(n + 1).join(' ') + '- ';
}

function getMapped(param: PD.Mapped<any>, offset: number) {
    let ret = '';
    for (const [n] of param.select.options) {
        ret += offsetS(offset);
        ret += `**${n}**:\n`;
        ret += paramInfo(param.map(n), offset + 2);
        ret += '\n';
    }
    return ret;
}

function getParams(params: PD.Params, offset: number) {
    let ret = '';
    for (const k of Object.keys(params)) {
        const param = params[k];
        ret += offsetS(offset);
        ret += `**${k}**${param.isOptional ? '?:' : ':'} ${paramInfo(param, offset)}`;
        // if (param.defaultValue) {
        //     ret += ` = ${JSON.stringify(param.defaultValue)}`;
        // }
        if (param.description) {
            ret += ` *(${param.description})*`;
        }
        ret += '\n';
    }
    return ret;
}