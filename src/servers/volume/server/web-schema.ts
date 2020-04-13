/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { VOLUME_SERVER_VERSION } from './version';
import { LimitsConfig, ServerConfig } from '../config';

export function getSchema() {
    function detail(i: number) {
        return `${i} (${Math.round(100 * LimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel[i] / 1000 / 1000) / 100 }M voxels)`;
    }
    const detailMax = LimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel.length - 1;
    const sources = ServerConfig.idMap.map(m => m[0]);

    return {
        openapi: '3.0.0',
        info: {
            version: VOLUME_SERVER_VERSION,
            title: 'Volume Server',
            description: 'The VolumeServer is a service for accessing subsets of volumetric data. It automatically downsamples the data depending on the volume of the requested region to reduce the bandwidth requirements and provide near-instant access to even the largest data sets.',
        },
        tags: [
            {
                name: 'General',
            }
        ],
        paths: {
            [`${ServerConfig.apiPrefix}/{source}/{id}/`]: {
                get: {
                    tags: ['General'],
                    summary: 'Returns a JSON response specifying if data is available and the maximum region that can be queried.',
                    operationId: 'getInfo',
                    parameters: [
                        { $ref: '#/components/parameters/source' },
                        { $ref: '#/components/parameters/id' },
                    ],
                    responses: {
                        200: {
                            description: 'Volume availability and info',
                            content: {
                                'application/json': {
                                    schema: { $ref: '#/components/schemas/info' }
                                }
                            }
                        },
                    },
                }
            },
            [`${ServerConfig.apiPrefix}/{source}/{id}/box/{a1,a2,a3}/{b1,b2,b3}/`]: {
                get: {
                    tags: ['General'],
                    summary: 'Returns density data inside the specified box for the given entry. For X-ray data, returns 2Fo-Fc and Fo-Fc volumes in a single response.',
                    operationId: 'getBox',
                    parameters: [
                        { $ref: '#/components/parameters/source' },
                        { $ref: '#/components/parameters/id' },
                        {
                            name: 'bottomLeftCorner',
                            in: 'path',
                            description: 'Bottom left corner of the query region in Cartesian or fractional coordinates (determined by the `space` query parameter).',
                            required: true,
                            schema: {
                                type: 'list',
                                items: {
                                    type: 'float',
                                }
                            },
                            style: 'simple'
                        },
                        {
                            name: 'topRightCorner',
                            in: 'path',
                            description: 'Top right corner of the query region in Cartesian or fractional coordinates (determined by the `space` query parameter).',
                            required: true,
                            schema: {
                                type: 'list',
                                items: {
                                    type: 'float',
                                }
                            },
                            style: 'simple'
                        },
                        { $ref: '#/components/parameters/encoding' },
                        { $ref: '#/components/parameters/detail' },
                        {
                            name: 'space',
                            in: 'query',
                            description: 'Determines the coordinate space the query is in. Can be cartesian or fractional. An optional argument, default values is cartesian.',
                            schema: {
                                type: 'string',
                                enum: ['cartesian', 'fractional']
                            },
                            style: 'form'
                        }
                    ],
                    responses: {
                        200: {
                            description: 'Volume box',
                            content: {
                                'text/plain': {},
                                'application/octet-stream': {},
                            }
                        },
                    },
                }
            },
            [`${ServerConfig.apiPrefix}/{source}/{id}/cell/`]: {
                get: {
                    tags: ['General'],
                    summary: 'Returns (downsampled) volume data for the entire "data cell". For X-ray data, returns unit cell of 2Fo-Fc and Fo-Fc volumes, for EM data returns everything.',
                    operationId: 'getCell',
                    parameters: [
                        { $ref: '#/components/parameters/source' },
                        { $ref: '#/components/parameters/id' },
                        { $ref: '#/components/parameters/encoding' },
                        { $ref: '#/components/parameters/detail' },
                    ],
                    responses: {
                        200: {
                            description: 'Volume cell',
                            content: {
                                'text/plain': {},
                                'application/octet-stream': {},
                            }
                        },
                    },
                }
            }
        },
        components: {
            schemas: {
                // TODO how to keep in sync with (or derive from) `api.ts/ExtendedHeader`
                info: {
                    properties: {
                        formatVersion: {
                            type: 'string',
                            description: 'Format version number'
                        },
                        axisOrder: {
                            type: 'array',
                            items: { type: 'number' },
                            description: 'Axis order from the slowest to fastest moving, same as in CCP4'
                        },
                        origin: {
                            type: 'array',
                            items: { type: 'number' },
                            description: 'Origin in fractional coordinates, in axisOrder'
                        },
                        dimensions: {
                            type: 'array',
                            items: { type: 'number' },
                            description: 'Dimensions in fractional coordinates, in axisOrder'
                        },
                        spacegroup: {
                            properties: {
                                number: { type: 'number' },
                                size: {
                                    type: 'array',
                                    items: { type: 'number' }
                                },
                                angles: {
                                    type: 'array',
                                    items: { type: 'number' }
                                },
                                isPeriodic: {
                                    type: 'boolean',
                                    description: 'Determine if the data should be treated as periodic or not. (e.g. X-ray = periodic, EM = not periodic)'
                                },
                            }
                        },
                        channels: {
                            type: 'array',
                            items: { type: 'string' }
                        },
                        valueType: {
                            type: 'string',
                            enum: ['float32', 'int16', 'int8'],
                            description: 'Determines the data type of the values'
                        },
                        blockSize: {
                            type: 'number',
                            description: 'The value are stored in blockSize^3 cubes'
                        },
                        sampling: {
                            type: 'array',
                            items: {
                                properties: {
                                    byteOffset: { type: 'number' },
                                    rate: {
                                        type: 'number',
                                        description: 'How many values along each axis were collapsed into 1'
                                    },
                                    valuesInfo: {
                                        properties: {
                                            mean: { type: 'number' },
                                            sigma: { type: 'number' },
                                            min: { type: 'number' },
                                            max: { type: 'number' },
                                        }
                                    },
                                    sampleCount: {
                                        type: 'array',
                                        items: { type: 'number' },
                                        description: 'Number of samples along each axis, in axisOrder'
                                    },
                                }
                            }
                        }
                    }
                }
            },
            parameters: {
                source: {
                    name: 'source',
                    in: 'path',
                    description: `Specifies the data source (determined by the experiment method). Currently supported sources are: ${sources.join(', ')}.`,
                    required: true,
                    schema: {
                        type: 'string',
                        enum: sources
                    },
                    style: 'simple'
                },
                id: {
                    name: 'id',
                    in: 'path',
                    description: 'Id of the entry. For x-ray, use PDB ID (i.e. 1cbs) and for em use EMDB id (i.e. emd-8116).',
                    required: true,
                    schema: {
                        type: 'string',
                    },
                    style: 'simple'
                },
                encoding: {
                    name: 'encoding',
                    in: 'query',
                    description: 'Determines if text based CIF or binary BinaryCIF encoding is used. An optional argument, default is BinaryCIF encoding.',
                    schema: {
                        type: 'string',
                        enum: ['cif', 'bcif']
                    },
                    style: 'form'
                },
                detail: {
                    name: 'detail',
                    in: 'query',
                    description: `Determines the maximum number of voxels the query can return. Possible values are in the range from ${detail(0)} to ${detail(detailMax)}. Default value is 0. Note: different detail levels might lead to the same result.`,
                    schema: {
                        type: 'integer',
                    },
                    style: 'form'
                }
            }
        }
    };
}

export const shortcutIconLink = `<link rel='shortcut icon' href='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAMAAABEpIrGAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAnUExURQAAAMIrHrspHr0oH7soILonHrwqH7onILsoHrsoH7soH7woILwpIKgVokoAAAAMdFJOUwAQHzNxWmBHS5XO6jdtAmoAAACZSURBVDjLxZNRCsQgDAVNXmwb9f7nXZEaLRgXloXOhwQdjMYYwpOLw55fBT46KhbOKhmRR2zLcFJQj8UR+HxFgArIF5BKJbEncC6NDEdI5SatBRSDJwGAoiFDONrEJXWYhGMIcRJGCrb1TOtDahfUuQXd10jkFYq0ViIrbUpNcVT6redeC1+b9tH2WLR93Sx2VCzkv/7NjfABxjQHksGB7lAAAAAASUVORK5CYII=' />`;