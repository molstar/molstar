/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from https://github.com/prisma/graphql-request, Copyright (c) 2017 Graphcool, MIT
 */

import { RuntimeContext } from '../mol-task';
import { AssetManager, Asset } from './assets';

type Variables = { [key: string]: any }

interface GraphQLError {
    message: string
    locations: { line: number, column: number }[]
    path: string[]
}

interface GraphQLResponse {
    data?: any
    errors?: GraphQLError[]
    extensions?: any
    status: number
    [key: string]: any
}

interface GraphQLRequestContext {
    query: string
    variables?: Variables
}

export class ClientError extends Error {
    response: GraphQLResponse
    request: GraphQLRequestContext

    constructor (response: GraphQLResponse, request: GraphQLRequestContext) {
        const message = `${ClientError.extractMessage(response)}: ${JSON.stringify({ response, request })}`;

        super(message);

        this.response = response;
        this.request = request;

        // this is needed as Safari doesn't support .captureStackTrace
        if (typeof Error.captureStackTrace === 'function') {
            Error.captureStackTrace(this, ClientError);
        }
    }

    private static extractMessage (response: GraphQLResponse): string {
        return response.errors ? response.errors[0].message : `GraphQL Error (Code: ${response.status})`;
    }
}

export class GraphQLClient {
    constructor(private url: string, private assetManager: AssetManager) { }

    async request(ctx: RuntimeContext, query: string, variables?: Variables): Promise<Asset.Wrapper<'json'>> {

        const body = JSON.stringify({ query, variables }, null, 2);
        const url = Asset.getUrlAsset(this.assetManager, this.url, body);
        const result = await this.assetManager.resolve(url, 'json').runInContext(ctx);

        if (!result.data.errors && result.data.data) {
            return {
                data: result.data.data,
                dispose: result.dispose
            };
        } else {
            const errorResult = typeof result.data === 'string' ? { error: result.data } : result.data;
            throw new ClientError({ ...errorResult }, { query, variables });
        }
    }
}