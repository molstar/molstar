
export interface TokenizerState {
    data: string

    position: number
    length: number

    currentLineNumber: number
    currentTokenStart: number
    currentTokenEnd: number

    currentTokenType?: number
}
