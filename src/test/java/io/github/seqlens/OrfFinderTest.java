package io.github.seqlens;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.Test;

class OrfFinderTest {

    // ATG AAA TTT TAG starting at position 3 (1-based)
    private static final String FORWARD = "CCATGAAATTTTAGCC";

    @Test
    void findsForwardStrandOrf() {
        List<Orf> orfs = OrfFinder.find(FORWARD, 0);
        Orf orf = orfs.stream().filter(o -> o.strand() == '+').findFirst().orElseThrow();
        assertEquals(3, orf.frame());
        assertEquals(3, orf.start());
        assertEquals(14, orf.end());
        assertEquals("MKF", orf.protein());
    }

    @Test
    void findsReverseStrandOrfWithForwardCoordinates() {
        String reverse = DnaTools.reverseComplement(FORWARD);
        List<Orf> orfs = OrfFinder.find(reverse, 0);
        Orf orf = orfs.stream().filter(o -> o.strand() == '-').findFirst().orElseThrow();
        assertEquals(3, orf.start());
        assertEquals(14, orf.end());
        assertEquals("MKF", orf.protein());
    }

    @Test
    void respectsMinimumLength() {
        assertTrue(OrfFinder.find(FORWARD, 4).isEmpty());
    }

    @Test
    void ignoresOrfWithoutStopCodon() {
        assertTrue(OrfFinder.find("ATGAAAAAAAAA", 0).isEmpty());
    }

    @Test
    void reportsLongestFirst() {
        String seq = "ATGTAA" + "ATGAAAAAATAA";
        List<Orf> orfs = OrfFinder.find(seq, 0);
        assertEquals("MKK", orfs.get(0).protein());
    }
}
