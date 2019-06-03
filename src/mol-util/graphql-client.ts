/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from https://github.com/prisma/graphql-request, Copyright (c) 2017 Graphcool, MIT
 */

import { RuntimeContext } from '../mol-task';

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
        const message = `${ClientError.extractMessage(response)}: ${JSON.stringify({ response, request })}`

        super(message)

        this.response = response
        this.request = request

        // this is needed as Safari doesn't support .captureStackTrace
        /* tslint:disable-next-line */
        if (typeof Error.captureStackTrace === 'function') {
            Error.captureStackTrace(this, ClientError)
        }
    }

    private static extractMessage (response: GraphQLResponse): string {
        try {
            return response.errors![0].message
        } catch (e) {
            return `GraphQL Error (Code: ${response.status})`
        }
    }
}

export class GraphQLClient {
    constructor(private url: string, private fetch: import('../mol-util/data-source').AjaxTask) {
        this.url = url
    }

    async request<T extends any>(ctx: RuntimeContext, query: string, variables?: Variables): Promise<T> {

        const body = JSON.stringify({
            query,
            variables: variables ? variables : undefined,
        })

        const resultStr = await this.fetch({ url: this.url, body }).runInContext(ctx)
        const result = JSON.parse(resultStr)

        if (!result.errors && result.data) {
            return result.data
        } else {
            const errorResult = typeof result === 'string' ? { error: result } : result
            throw new ClientError(
                { ...errorResult },
                { query, variables },
            )
        }
    }
}