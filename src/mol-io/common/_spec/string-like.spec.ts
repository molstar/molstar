import { ChunkedBigString } from '../string-like';


const SAMPLE_ASCII = `Basketball is a team sport in which two teams, most commonly 
of five players each, opposing one another on a rectangular court, compete with 
the primary objective of shooting a basketball (approximately 9.4 inches (24 cm) 
in diameter) through the defender's hoop (a basket 18 inches (46 cm) in diameter 
mounted 10 feet (3.048 m) high to a backboard at each end of the court), while 
preventing the opposing team from shooting through their own hoop. A field 
goal is worth two points, unless made from behind the three-point line, when 
it is worth three. After a foul, timed play stops and the player fouled or 
designated to shoot a technical foul is given one, two or three one-point free 
throws. The team with the most points at the end of the game wins, but if 
regulation play expires with the score tied, an additional period of play 
(overtime) is mandated.
`;

const SAMPLE_UNICODE_CS = `Zámecký pivovar Frýdlant je průmyslový podnik ve Frýdlantě, 
městě na severu České republiky ve Frýdlantském výběžku Libereckého kraje. Pivovar 
se nachází pod místním hradem a zámkem na břehu řeky Smědé. Stojí na počátku staré 
cesty do Raspenavy (dnes silnice II/290). První zmínky o pivovaru jsou z konce 14. 
století, kdy se v něm vařilo pivo pro potřeby zdejšího hradu. O dvě století 
později se ovšem kvůli lepší dostupnosti přírodních zdrojů přestěhoval do 
podhradí na břeh řeky Smědé. V dalších letech se pivovar postupně rozšiřoval 
a zvyšovala se také jeho roční produkce, která na počátku 20. století 
přesahovala 30 tisíc hektolitrů piva. Ovšem s ohledem na objem finančních 
prostředků, jež by si vyžádala jeho modernizace, došlo v roce 1949 k ukončení 
provozu pivovaru a jeho prostory se využívaly pro zrání sýrů a uskladnění 
zeleniny. Na počátku 21. století se však objevily snahy o opětovné zahájení 
provozu. Vaří se zde pivo značky Albrecht. Slavnostní otevření pivovaru spolu 
s jeho požehnáním proběhlo v červenci 2014.
`;

const SAMPLE_UNICODE_PA = `ਵਿਕੀਪੀਡੀਆ, ਇੱਕ ਆਜ਼ਾਦ ਵਿਸ਼ਵਕੋਸ਼ ਤੋਂ 
ਬਾਸਕਟਬਾਲ 
ਮਾਈਕਲ ਜਾਰਡਨ ਪੁਰਾਣੇ ਬੌਸਟਨ ਗਾਰਡਨ ਵਿਖੇ ਸਲੈਮ ਡੰਕ ਮਾਰਨ ਜਾਂਦਾ ਹੋਇਆ
ਸਰਬ-ਉੱਚ ਅਦਾਰਾ	ਫ਼ੀਬਾ
ਪਹਿਲੋਂ ਖੇਡੀ ਗਈ	1891, ਸਪਰਿੰਗਫ਼ੀਲਡ, ਮੈਸਾਚੂਸਟਸ, ਅਮਰੀਕਾ
ਗੁਣ
ਛੋਹ	ਮੌਜੂਦ
ਜੁੱਟ ਵਿੱਚ ਜੀਅ	ਇੱਕ ਪਾਸੇ 5
ਰਲ਼ਵਾਂ ਲਿੰਗ	ਹਾਂ, ਅੱਡੋ-ਅੱਡ ਮੁਕਾਬਲੇ
ਕਿਸਮ	ਜੁੱਟ ਖੇਡ, ਖਿੱਦੋ ਖੇਡ
ਸਾਜ਼ੋ-ਸਮਾਨ	ਬਾਸਕਟਬਾਲ
ਟਿਕਾਣਾ	ਅੰਦਰਲੇ ਮੈਦਾਨ (ਮੁੱਖ) ਜਾਂ ਬਾਹਰਲੇ ਮੈਦਾਨ (ਸਟਰੀਟਬਾਲ)
ਮੌਜੂਦਗੀ
ਓਲੰਪਿਕ	1904 ਅਤੇ 1924 ਦੀਆਂ ਉਲੰਪਿਕ ਖੇਡਾਂ 'ਚ ਵਿਖਾਈ ਗਈ
1936 ਤੋਂ ਗਰਮੀਆਂ ਦੀ ਓਲੰਪਿਕ ਦਾ ਹਿੱਸਾ
ਬਾਸਕਟਬਾਲ ਪੰਜ ਖਿਡਾਰੀਆਂ ਦੇ ਦੋ ਜੁੱਟਾਂ ਵੱਲੋਂ ਕਿਸੇ ਚੌਭੁਜੀ ਮੈਦਾਨ ਉੱਤੇ ਖੇਡੀ ਜਾਣ ਵਾਲ਼ੀ 
ਇੱਕ ਖੇਡ ਹੈ। ਮੁੱਖ ਮਕਸਦ ਦੋਹੇਂ ਸਿਰਿਆਂ ਉੱਤੇ ਗੱਡੇ ਇੱਕ ਖੰਭੇ ਉੱਤੇ ਲੱਗੀ 10 ਫੁੱਟ (3 
ਮੀ.) ਉੱਚੀ ਅਤੇ 18 ਇੰਚ (46 ਸੈ.ਮੀ.) ਦੇ ਵਿਆਸ ਵਾਲ਼ੀ ਬਿਨਾਂ ਤਲੇ ਵਾਲ਼ੀ 
ਜਾਲ਼ੀਦਾਰ ਟੋਕਰੀ ਵਿੱਚ ਗੇਂਦ ਮਾਰਨਾ ਹੁੰਦਾ ਹੈ।ਬਾਸਕਟਬਾਲ ਦੇ ਗ੍ਰਾਉੰਡ ਦੀ ਲੰਬਾਈ 28 
ਮੀ: ਤੇ ਚੌੜਾਈ 15 ਮੀ: ਹੁੰਦੀ ਹੈ।ਬਾਸਕਟਬਾਲ ਦੁਨੀਆ ਦੀਆਂ ਸਭ ਤੋਂ ਮਸ਼ਹੂਰ ਅਤੇ 
ਮਕਬੂਲ ਖੇਡਾਂ ਵਿੱਚੋਂ ਇੱਕ ਹੈ।[1]

ਬਾਸਕਟਬਾਲ ਇੱਕ ਟੀਮ ਖੇਡ ਹੈ ਜਿਸ ਵਿੱਚ ਦੋ ਟੀਮਾਂ, ਆਮ ਤੌਰ 'ਤੇ ਪੰਜ ਖਿਡਾਰੀ, 
ਇੱਕ ਆਇਤਾਕਾਰ ਅਦਾਲਤ ਵਿੱਚ ਇੱਕ ਦੂਜੇ ਦਾ ਵਿਰੋਧ ਕਰਨ ਵਾਲੀਆਂ, ਇੱਕ 
ਬਾਸਕਟਬਾਲ (ਲਗਭਗ 9.4 ਇੰਚ (24 ਸੈ) ਵਿਆਸ) ਵਿੱਚ ਨਿਸ਼ਾਨਾ ਲਗਾਉਣ ਦੇ 
ਮੁੱਢਲੇ ਉਦੇਸ਼ ਨਾਲ ਮੁਕਾਬਲਾ ਕਰਦੇ ਹਨ। ਇੱਕ ਟੋਕਰੀ 18 ਇੰਚ (46 ਸੈਂਟੀਮੀਟਰ) 
ਵਿਆਸ ਵਾਲੀ ਇੱਕ ਫੁੱਟ 10 ਫੁੱਟ (3.048 ਮੀਟਰ) ਉੱਚੀ ਇੱਕ ਅਦਾਲਤ ਦੇ ਹਰ 
ਸਿਰੇ 'ਤੇ ਇੱਕ ਬਕਬੋਰਡ ਤੇ ਪਈ) ਜਦੋਂ ਕਿ ਵਿਰੋਧੀ ਟੀਮ ਨੂੰ ਉਨ੍ਹਾਂ ਦੇ ਆਪਣੇ ਹੂਪ 
ਦੁਆਰਾ ਗੋਲੀ ਮਾਰਨ ਤੋਂ ਰੋਕਿਆ। ਇੱਕ ਫੀਲਡ ਟੀਚਾ ਦੋ ਪੁਆਇੰਟਾਂ ਦਾ ਮੁੱਲਵਾਨ ਹੁੰਦਾ 
ਹੈ, ਜਦੋਂ ਤੱਕ ਕਿ ਤਿੰਨ-ਪੁਆਇੰਟ ਦੀ ਰੇਖਾ ਦੇ ਪਿੱਛੇ ਨਹੀਂ ਬਣਾਇਆ ਜਾਂਦਾ, ਜਦੋਂ ਇਹ 
ਤਿੰਨ ਦੀ ਕੀਮਤ ਵਾਲਾ ਹੁੰਦਾ ਹੈ। ਇੱਕ ਅਸ਼ੁੱਧ ਦੇ ਬਾਅਦ, ਸਮੇਂ ਸਿਰ ਖੇਡ ਰੁਕ ਜਾਂਦੀ ਹੈ 
ਅਤੇ ਖਿਡਾਰੀ ਨੂੰ ਤਕਨੀਕੀ ਫਾ .ਲ ਸ਼ੂਟ ਕਰਨ ਲਈ ਨਿਰਧਾਰਤ ਕੀਤਾ ਜਾਂਦਾ ਹੈ ਜਿਸ ਨੂੰ 
ਇੱਕ ਜਾਂ ਵਧੇਰੇ ਇੱਕ-ਪੁਆਇੰਟ ਮੁਫਤ ਥ੍ਰੋਅ ਦਿੱਤਾ ਜਾਂਦਾ ਹੈ। ਖੇਡ ਦੇ ਅੰਤ ਵਿੱਚ ਸਭ ਤੋਂ 
ਜ਼ਿਆਦਾ ਪੁਆਇੰਟਾਂ ਵਾਲੀ ਟੀਮ ਜਿੱਤ ਜਾਂਦੀ ਹੈ, ਪਰ ਜੇ ਨਿਯਮਿਤ ਖੇਡ ਸਕੋਰ ਦੇ ਬਰਾਬਰੀ 
ਨਾਲ ਖਤਮ ਹੋ ਜਾਂਦੀ ਹੈ, ਤਾਂ ਵਾਧੂ ਸਮੇਂ ਦਾ ਖੇਡ (ਓਵਰਟਾਈਮ) ਲਾਜ਼ਮੀ ਹੁੰਦਾ ਹੈ।
`;


const TESTING_LOG_STRING_CHUNK_SIZE = 3; // chunk of size 8

function testUtf8Decoding(text: string) {
    const bytes = Buffer.from(text, 'utf-8');
    const bigString = ChunkedBigString.fromUtf8Data(bytes, undefined, undefined, TESTING_LOG_STRING_CHUNK_SIZE);
    const redecoded = bigString.toString();
    expect(redecoded).toEqual(text);
}

describe('ChunkedBigString.fromUtf8Data', () => {
    test('decode ASCII', async () => {
        testUtf8Decoding(SAMPLE_ASCII);
    });

    test('decode CS', async () => {
        testUtf8Decoding(SAMPLE_UNICODE_CS);
    });

    test('decode PA', async () => {
        testUtf8Decoding(SAMPLE_UNICODE_PA);
    });
});

describe('ChunkedBigString.at', () => {
    test('at ASCII', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_ASCII, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.at(0))
            .toEqual('B');
        expect(bigString.at(SAMPLE_ASCII.indexOf('9')))
            .toEqual('9');
        expect(bigString.at(10_000))
            .toEqual(undefined);
    });

    test('at CS', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.at(0))
            .toEqual('Z');
        expect(bigString.at(SAMPLE_UNICODE_CS.indexOf('ř')))
            .toEqual('ř');
        expect(bigString.at(10_000))
            .toEqual(undefined);
        expect(bigString.at(-10_000))
            .toEqual(undefined);
        expect(bigString.at(-10))
            .toEqual('n');
    });
});

describe('ChunkedBigString.charAt', () => {
    test('charAt ASCII', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_ASCII, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.charAt(0))
            .toEqual('B');
        expect(bigString.charAt(SAMPLE_ASCII.indexOf('9')))
            .toEqual('9');
        expect(bigString.charAt(10_000))
            .toEqual('');
    });

    test('charAt CS', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.charAt(0))
            .toEqual('Z');
        expect(bigString.charAt(SAMPLE_UNICODE_CS.indexOf('ř')))
            .toEqual('ř');
        expect(bigString.charAt(10_000))
            .toEqual('');
        expect(bigString.charAt(-10_000))
            .toEqual('');
        expect(bigString.charAt(-10))
            .toEqual('');
    });
});

describe('ChunkedBigString.charCodeAt', () => {
    test('charCodeAt ASCII', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_ASCII, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.charCodeAt(0))
            .toEqual('B'.charCodeAt(0));
        expect(bigString.charCodeAt(SAMPLE_ASCII.indexOf('9')))
            .toEqual('9'.charCodeAt(0));
        expect(bigString.charCodeAt(10_000))
            .toEqual(NaN);
    });

    test('charCodeAt CS', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.charCodeAt(0))
            .toEqual('Z'.charCodeAt(0));
        expect(bigString.charCodeAt(SAMPLE_UNICODE_CS.indexOf('ř')))
            .toEqual('ř'.charCodeAt(0));
        expect(bigString.charCodeAt(10_000))
            .toEqual(NaN);
    });
});

describe('ChunkedBigString.indexOf', () => {
    test('indexOf ASCII, within chunk', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_ASCII, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.indexOf('Basket'))
            .toEqual(SAMPLE_ASCII.indexOf('Basket'));
        expect(bigString.indexOf('another'))
            .toEqual(SAMPLE_ASCII.indexOf('another'));
        expect(bigString.indexOf('unicorn'))
            .toEqual(-1);
    });

    test('indexOf ASCII, across chunks', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_ASCII, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.indexOf('sport'))
            .toEqual(SAMPLE_ASCII.indexOf('sport'));
    });

    test('indexOf CS, within chunk', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.indexOf('Zámecký'))
            .toEqual(SAMPLE_UNICODE_CS.indexOf('Zámecký'));
        expect(bigString.indexOf('zámkem'))
            .toEqual(SAMPLE_UNICODE_CS.indexOf('zámkem'));
        expect(bigString.indexOf('unicorn'))
            .toEqual(-1);
    });

    test('indexOf CS, across chunks', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.indexOf('městě'))
            .toEqual(SAMPLE_UNICODE_CS.indexOf('městě'));
    });

    test('indexOf CS, with position param', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.indexOf('pivovar')).toEqual(8);
        expect(bigString.indexOf('pivovar', 8)).toEqual(8);
        expect(bigString.indexOf('pivovar', 9)).toEqual(286);
        expect(bigString.indexOf('pivovar', 286)).toEqual(286);
        expect(bigString.indexOf('pivovar', 287)).toEqual(514);
        expect(bigString.indexOf('pivovar', 514)).toEqual(514);
        expect(bigString.indexOf('pivovar', 515)).toEqual(776);
        expect(bigString.indexOf('pivovar', 776)).toEqual(776);
        expect(bigString.indexOf('pivovar', 777)).toEqual(983);
        expect(bigString.indexOf('pivovar', 983)).toEqual(983);
        expect(bigString.indexOf('pivovar', 984)).toEqual(-1);
        expect(bigString.indexOf('pivovar', 10_000)).toEqual(-1);
        expect(bigString.indexOf('pivovar', -10_000)).toEqual(8);
    });
});

describe('ChunkedBigString.substring', () => {
    test('substring ASCII', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_ASCII, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.substring())
            .toEqual(SAMPLE_ASCII);
        expect(bigString.substring(0, 10))
            .toEqual('Basketball');
        expect(bigString.substring(SAMPLE_ASCII.indexOf('mandated')))
            .toEqual('mandated.\n');
        expect(bigString.substring(SAMPLE_ASCII.indexOf('opposing team'), SAMPLE_ASCII.indexOf('opposing team') + 'opposing team'.length))
            .toEqual('opposing team');
    });

    test('substring CS', async () => {
        const bigString = ChunkedBigString.fromString(SAMPLE_UNICODE_CS, TESTING_LOG_STRING_CHUNK_SIZE);
        expect(bigString.substring())
            .toEqual(SAMPLE_UNICODE_CS);
        expect(bigString.substring(0, 7))
            .toEqual('Zámecký');
        expect(bigString.substring(SAMPLE_UNICODE_CS.indexOf('2014')))
            .toEqual('2014.\n');
        expect(bigString.substring(SAMPLE_UNICODE_CS.indexOf('Slavnostní otevření'), SAMPLE_UNICODE_CS.indexOf('Slavnostní otevření') + 'Slavnostní otevření'.length))
            .toEqual('Slavnostní otevření');
    });
});

// TODO enforce short chunks (2**3) for tests, otherwise they don't make sense
