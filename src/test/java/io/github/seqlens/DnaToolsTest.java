package io.github.seqlens;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Map;
import org.junit.jupiter.api.Test;

class DnaToolsTest {

    @Test
    void gcContentOfPureGcIsOne() {
        assertEquals(1.0, DnaTools.gcContent("GGCC"), 1e-9);
    }

    @Test
    void gcContentIgnoresAmbiguousBases() {
        assertEquals(0.5, DnaTools.gcContent("ATGCNNNN"), 1e-9);
    }

    @Test
    void gcContentOfEmptySequenceIsZero() {
        assertEquals(0.0, DnaTools.gcContent(""), 1e-9);
    }

    @Test
    void reverseComplement() {
        assertEquals("NGCAT", DnaTools.reverseComplement("ATGCN"));
        assertEquals("ACGT", DnaTools.reverseComplement("ACGT")); // palindrome
    }

    @Test
    void reverseComplementTwiceGivesOriginal() {
        String seq = "ATGCGTACGTTAGCRYKM";
        assertEquals(seq, DnaTools.reverseComplement(DnaTools.reverseComplement(seq)));
    }

    @Test
    void transcribeReplacesThymine() {
        assertEquals("AUGCUU", DnaTools.transcribe("atgctt"));
    }

    @Test
    void baseCounts() {
        Map<Character, Integer> c = DnaTools.baseCounts("AACGTTTN");
        assertEquals(2, c.get('A'));
        assertEquals(1, c.get('C'));
        assertEquals(1, c.get('G'));
        assertEquals(3, c.get('T'));
        assertEquals(1, c.get('N'));
    }

    @Test
    void meltingTemperatureWallaceRuleForShortOligos() {
        assertEquals(12.0, DnaTools.meltingTemperature("ATGC"), 1e-9);
    }

    @Test
    void kmerCountsAreOverlapping() {
        Map<String, Integer> k = DnaTools.kmerCounts("AAAA", 2);
        assertEquals(Map.of("AA", 3), k);
    }

    @Test
    void kmerCountsSkipAmbiguousKmers() {
        Map<String, Integer> k = DnaTools.kmerCounts("ACNGT", 2);
        assertEquals(Map.of("AC", 1, "GT", 1), k);
    }

    @Test
    void kmerRejectsNonPositiveK() {
        assertThrows(IllegalArgumentException.class, () -> DnaTools.kmerCounts("ACGT", 0));
    }

    @Test
    void validDna() {
        assertTrue(DnaTools.isValidDna("acgtNRY"));
        assertFalse(DnaTools.isValidDna("ACGTXZ"));
    }
}
