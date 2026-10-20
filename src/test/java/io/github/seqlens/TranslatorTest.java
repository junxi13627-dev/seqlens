package io.github.seqlens;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

class TranslatorTest {

    @Test
    void translatesStandardCode() {
        assertEquals("MAIVMGR*", Translator.translate("ATGGCCATTGTAATGGGCCGCTGA"));
    }

    @Test
    void translatesInOtherFrames() {
        // frame 2 skips one base, frame 3 skips two
        assertEquals("M", Translator.translate("CATG", 1));
        assertEquals("M", Translator.translate("CCATG", 2));
    }

    @Test
    void acceptsRnaAndLowercase() {
        assertEquals('M', Translator.translateCodon("aug"));
        assertEquals('*', Translator.translateCodon("UAA"));
    }

    @Test
    void allStopCodons() {
        for (String stop : new String[] {"TAA", "TAG", "TGA"}) {
            assertEquals('*', Translator.translateCodon(stop));
        }
    }

    @Test
    void ambiguousCodonBecomesX() {
        assertEquals('X', Translator.translateCodon("ANG"));
    }

    @Test
    void rejectsInvalidFrame() {
        assertThrows(IllegalArgumentException.class, () -> Translator.translate("ATG", 3));
    }
}
