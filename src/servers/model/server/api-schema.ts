/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import VERSION from '../version';
import { QueryParamInfo, QueryParamType, QueryDefinition, CommonQueryParamsInfo, QueryList } from './api';
import { ModelServerConfig as ServerConfig } from '../config';
import { MultipleQuerySpec } from './api-web-multiple';

export const shortcutIconLink = `<link rel='shortcut icon' href='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAMAAABEpIrGAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAnUExURQAAAMIrHrspHr0oH7soILonHrwqH7onILsoHrsoH7soH7woILwpIKgVokoAAAAMdFJOUwAQHzNxWmBHS5XO6jdtAmoAAACZSURBVDjLxZNRCsQgDAVNXmwb9f7nXZEaLRgXloXOhwQdjMYYwpOLw55fBT46KhbOKhmRR2zLcFJQj8UR+HxFgArIF5BKJbEncC6NDEdI5SatBRSDJwGAoiFDONrEJXWYhGMIcRJGCrb1TOtDahfUuQXd10jkFYq0ViIrbUpNcVT6redeC1+b9tH2WLR93Sx2VCzkv/7NjfABxjQHksGB7lAAAAAASUVORK5CYII=' />`;

export function getApiSchema() {
    return {
        openapi: '3.0.0',
        info: {
            version: VERSION,
            title: 'ModelServer',
            description: 'The ModelServer is a service for accessing subsets of macromolecular model data.',
        },
        tags: [
            {
                name: 'General',
            }
        ],
        paths: getPaths(),
        components: {
            parameters: {
                id: {
                    name: 'id',
                    in: 'path',
                    description: 'Id of the entry (i.e. 1tqn).',
                    required: true,
                    schema: {
                        type: 'string',
                    },
                    style: 'simple'
                },
            }
        }
    };
}

function getPaths() {
    const ret: any = {};
    for (const { name, definition } of QueryList) {
        ret[`${ServerConfig.apiPrefix}/v1/{id}/${name}`] = getQueryInfo(definition);
    }

    const queryManySummary = 'Executes multiple queries at the same time and writes them as separate data blocks.';
    const queryManyExample: MultipleQuerySpec = {
        queries: [
            { entryId: '1cbs', query: 'residueInteraction', params: { atom_site: [{ label_comp_id: 'REA' }], radius: 5 } },
            { entryId: '1tqn', query: 'full', copy_all_categories: true }
        ],
        encoding: 'cif',
        asTarGz: false
    };
    ret[`${ServerConfig.apiPrefix}/v1/query-many`] = {
        get: {
            tags: ['General'],
            summary: queryManySummary,
            operationId: 'query-many',
            parameters: [{
                name: 'query',
                in: 'query',
                description: 'URI encoded JSON object with the query definiton.',
                required: true,
                schema: {
                    type: 'string',
                },
                example: JSON.stringify(queryManyExample),
                style: 'simple'
            }],
            responses: {
                200: {
                    description: 'Separate CIF data blocks with the result',
                    content: {
                        'text/plain': {},
                        'application/octet-stream': {},
                    }
                }
            }
        },
        post: {
            tags: ['General'],
            summary: queryManySummary,
            operationId: 'query-many-post',
            parameters: [],
            requestBody: {
                content: {
                    'application/json': {
                        schema: { type: 'object' },
                        example: queryManyExample
                    }
                }
            },
            responses: {
                200: {
                    description: 'Separate CIF data blocks with the result',
                    content: {
                        'text/plain': {},
                        'application/octet-stream': {},
                    }
                }
            }
        }
    };

    return ret;
}

function getQueryInfo(def: QueryDefinition) {
    const jsonExample: any = {};
    def.jsonParams.forEach(p => {
        if (!p.exampleValues || !p.exampleValues.length) return;
        jsonExample[p.name] = p.exampleValues[0];
    });

    return {
        get: {
            tags: ['General'],
            summary: def.description,
            operationId: def.name,
            parameters: [
                { $ref: '#/components/parameters/id' },
                ...def.restParams.map(getParamInfo),
                ...CommonQueryParamsInfo.map(getParamInfo)
            ],
            responses: {
                200: {
                    description: def.description,
                    content: {
                        'text/plain': {},
                        'application/octet-stream': {},
                    }
                }
            }
        },
        post: {
            tags: ['General'],
            summary: def.description,
            operationId: def.name + '-post',
            parameters: [
                { $ref: '#/components/parameters/id' },
                ...CommonQueryParamsInfo.map(getParamInfo)
            ],
            requestBody: {
                content: {
                    'application/json': {
                        schema: { type: 'object' },
                        example: jsonExample
                    }
                }
            },
            responses: {
                200: {
                    description: def.description,
                    content: {
                        'text/plain': {},
                        'application/octet-stream': {},
                    }
                }
            }
        }
    };
}

function getParamInfo(info: QueryParamInfo) {
    return {
        name: info.name,
        in: 'query',
        description: info.description,
        required: !!info.required,
        schema: {
            type: info.type === QueryParamType.String
                ? 'string' : info.type === QueryParamType.Integer
                    ? 'integer'
                    : info.type === QueryParamType.Boolean
                        ? 'boolean'
                        : 'number',
            enum: info.supportedValues ? info.supportedValues : void 0,
            default: info.defaultValue
        },
        style: 'simple'
    };
}