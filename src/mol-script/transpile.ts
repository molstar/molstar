import Transpiler from './transpilers/transpiler'
import _transpiler from './transpilers/all'
const transpiler: {[index: string]: Transpiler} = _transpiler

export function parse(lang: string, str: string) {
    try {
	const query = transpiler[lang](str);
	return query;      
    } catch (e) {

    }
}
