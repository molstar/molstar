/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import VERSION from './version';
import { Config } from './config';

export function getSchema(config: Config) {
    function mapPath(path: string) {
        return `${config.api_prefix}/${path}`;
    }

    return {
        openapi: '3.0.0',
        info: {
            version: VERSION,
            title: 'PluginState Server',
            description: 'The PluginState Server is a simple service for storing and retreiving states of the Mol* Viewer app.',
        },
        tags: [
            {
                name: 'General',
            }
        ],
        paths: {
            [mapPath(`list/`)]: {
                get: {
                    tags: ['General'],
                    summary: 'Returns a JSON response with the list of currently stored states.',
                    operationId: 'list',
                    parameters: [],
                    responses: {
                        200: {
                            description: 'A list of stored states',
                            content: {
                                'application/json': { }
                            }
                        }
                    },
                }
            },
            [mapPath(`get/{id}`)]: {
                get: {
                    tags: ['General'],
                    summary: 'Returns the Mol* Viewer state with the given id.',
                    operationId: 'get',
                    parameters: [
                        {
                            name: 'id',
                            in: 'path',
                            description: `Id of the state.`,
                            required: true,
                            schema: { type: 'string' },
                            style: 'simple'
                        }
                    ],
                    responses: {
                        200: {
                            description: 'A JSON object with the state.',
                            content: {
                                'application/json': { }
                            }
                        }
                    },
                }
            },
            [mapPath(`remove/{id}`)]: {
                get: {
                    tags: ['General'],
                    summary: 'Removes the Mol* Viewer state with the given id.',
                    operationId: 'remove',
                    parameters: [
                        {
                            name: 'id',
                            in: 'path',
                            description: `Id of the state.`,
                            required: true,
                            schema: { type: 'string' },
                            style: 'simple'
                        }
                    ],
                    responses: {
                        200: {
                            description: 'Empty response.',
                            content: { 'text/plain': { } }
                        }
                    },
                }
            },
            [mapPath(`set/`)]: {
                post: {
                    tags: ['General'],
                    summary: `Post Mol* Viewer state to the server. At most ${config.max_states} states can be stored. If the limit is reached, older states will be removed.`,
                    operationId: 'set',
                    requestBody: {
                        content: {
                            'application/json': {
                                schema: { type: 'object' }
                            }
                        }
                    },
                    parameters: [
                        {
                            name: 'name',
                            in: 'query',
                            description: `Name of the state. If none provided, current UTC date-time is used.`,
                            required: false,
                            schema: { type: 'string' },
                            style: 'simple'
                        },
                        {
                            name: 'description',
                            in: 'query',
                            description: `Description of the state.`,
                            required: false,
                            schema: { type: 'string' },
                            style: 'simple'
                        }
                    ],
                    responses: {
                        200: {
                            description: 'Empty response.',
                            content: { 'text/plain': { } }
                        }
                    },
                }
            },
        }
    };
}

export const shortcutIconLink = `<link rel='shortcut icon' href='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAMAAABEpIrGAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAnUExURQAAAMIrHrspHr0oH7soILonHrwqH7onILsoHrsoH7soH7woILwpIKgVokoAAAAMdFJOUwAQHzNxWmBHS5XO6jdtAmoAAACZSURBVDjLxZNRCsQgDAVNXmwb9f7nXZEaLRgXloXOhwQdjMYYwpOLw55fBT46KhbOKhmRR2zLcFJQj8UR+HxFgArIF5BKJbEncC6NDEdI5SatBRSDJwGAoiFDONrEJXWYhGMIcRJGCrb1TOtDahfUuQXd10jkFYq0ViIrbUpNcVT6redeC1+b9tH2WLR93Sx2VCzkv/7NjfABxjQHksGB7lAAAAAASUVORK5CYII=' />`;