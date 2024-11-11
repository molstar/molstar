/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { VERSION } from './version';
import { MembraneServerConfig } from './config';

export function getSchema() {
    return {
        openapi: '3.0.0',
        info: {
            version: VERSION,
            title: 'Membrane Server',
            description: 'The MembraneServer process an entry and predicts the orientation of the membrane layer, which can be used to compose molecular scenes using MolViewSpec.',
        },
        tags: [
            {
                name: 'General',
            }
        ],
        paths: {
            [`${MembraneServerConfig.apiPrefix}/predict/{id}/`]: {
                get: {
                    tags: ['General'],
                    summary: 'Returns a JSON response specifying if data is available and the maximum region that can be queried.',
                    operationId: 'predictMembrane',
                    parameters: [
                        { $ref: '#/components/parameters/id' },
                        { $ref: '#/components/parameters/assemblyId' },
                        { $ref: '#/components/parameters/numberOfSpherePoints' },
                        { $ref: '#/components/parameters/stepSize' },
                        { $ref: '#/components/parameters/minThickness' },
                        { $ref: '#/components/parameters/maxThickness' },
                        { $ref: '#/components/parameters/asaCutoff' },
                        { $ref: '#/components/parameters/adjust' },
                        { $ref: '#/components/parameters/tmdetDefinition' },
                    ],
                    responses: {
                        200: {
                            description: '',
                            content: {
                                'application/json': {
                                    schema: { $ref: '#/components/schemas/prediction' }
                                }
                            }
                        },
                    },
                }
            },
        },
        components: {
            schemas: {
                prediction: {
                    type: 'object',
                    properties: {
                        planePoint1: {
                            type: 'array',
                            items: {
                                type: 'number'
                            },
                            minItems: 3,
                            maxItems: 3,
                            description: 'Array of three numbers representing the first plane point'
                        },
                        planePoint2: {
                            type: 'array',
                            items: {
                                type: 'number'
                            },
                            minItems: 3,
                            maxItems: 3,
                            description: 'Array of three numbers representing the second plane point'
                        },
                        normalVector: {
                            type: 'array',
                            items: {
                                type: 'number'
                            },
                            minItems: 3,
                            maxItems: 3,
                            description: 'Array of three numbers representing the normal vector'
                        },
                        centroid: {
                            type: 'array',
                            items: {
                                type: 'number'
                            },
                            minItems: 3,
                            maxItems: 3,
                            description: 'Array of three numbers representing the centroid'
                        },
                        radius: {
                            type: 'number',
                            description: 'A number representing the radius'
                        }
                    }
                },
            },
            parameters: {
                id: {
                    name: 'id',
                    in: 'path',
                    description: 'Entry identifier of the entry, e.g. 1brr.',
                    required: true,
                    schema: {
                        default: '5cbg',
                        type: 'string',
                    },
                    style: 'simple',
                },
                assemblyId: {
                    name: 'assemblyId',
                    in: 'query',
                    description: 'Assembly identifier, e.g. 1',
                    required: false,
                    schema: {
                        default: '1',
                        type: 'string',
                    },
                    style: 'simple'
                },
                numberOfSpherePoints: {
                    name: 'numberOfSpherePoints',
                    in: 'query',
                    description: 'Number of spheres/directions to test for membrane placement. Original value is 350.',
                    required: false,
                    schema: {
                        type: 'integer',
                        minimum: 35,
                        maximum: 700,
                        default: 175,
                    },
                    style: 'simple'
                },
                stepSize: {
                    name: 'stepSize',
                    in: 'query',
                    description: 'Thickness of membrane slices that will be tested',
                    required: false,
                    schema: {
                        type: 'number',
                        minimum: 0.25,
                        maximum: 4,
                        default: 1,
                    },
                    style: 'simple'
                },
                minThickness: {
                    name: 'minThickness',
                    in: 'query',
                    description: 'Minimum membrane thickness used during refinement',
                    required: false,
                    schema: {
                        type: 'number',
                        minimum: 10,
                        maximum: 30,
                        default: 20,
                    },
                    style: 'simple'
                },
                maxThickness: {
                    name: 'maxThickness',
                    in: 'query',
                    description: 'Maximum membrane thickness used during refinement',
                    required: false,
                    schema: {
                        type: 'integer',
                        minimum: 30,
                        maximum: 50,
                        default: 40,
                    },
                    style: 'simple'
                },
                asaCutoff: {
                    name: 'asaCutoff',
                    in: 'query',
                    description: 'Relative ASA cutoff above which residues will be considered',
                    required: false,
                    schema: {
                        type: 'number',
                        minimum: 10,
                        maximum: 100,
                        default: 40,
                    },
                    style: 'simple'
                },
                adjust: {
                    name: 'adjust',
                    in: 'query',
                    description: 'Minimum length of membrane-spanning regions (original values: 14 for alpha-helices and 5 for beta sheets). Set to 0 to not optimize membrane thickness.',
                    required: false,
                    schema: {
                        type: 'integer',
                        minimum: 0,
                        maximum: 30,
                        default: 14,
                    },
                    style: 'simple'
                },
                tmdetDefinition: {
                    name: 'tmdetDefinition',
                    in: 'query',
                    description: `Use TMDET's classification of membrane-favoring amino acids. TMDET's classification shows better performance on porins and other beta-barrel structures.`,
                    required: false,
                    schema: {
                        type: 'boolean',
                        default: false,
                    },
                    style: 'simple'
                },
            }
        }
    };
}

export const shortcutIconLink = `<link rel='shortcut icon' href='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAMAAABEpIrGAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAnUExURQAAAMIrHrspHr0oH7soILonHrwqH7onILsoHrsoH7soH7woILwpIKgVokoAAAAMdFJOUwAQHzNxWmBHS5XO6jdtAmoAAACZSURBVDjLxZNRCsQgDAVNXmwb9f7nXZEaLRgXloXOhwQdjMYYwpOLw55fBT46KhbOKhmRR2zLcFJQj8UR+HxFgArIF5BKJbEncC6NDEdI5SatBRSDJwGAoiFDONrEJXWYhGMIcRJGCrb1TOtDahfUuQXd10jkFYq0ViIrbUpNcVT6redeC1+b9tH2WLR93Sx2VCzkv/7NjfABxjQHksGB7lAAAAAASUVORK5CYII=' />`;