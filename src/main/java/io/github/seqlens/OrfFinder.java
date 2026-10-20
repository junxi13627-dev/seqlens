package io.github.seqlens;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Finds open reading frames (ATG ... stop) in all six reading frames.
 *
 * <p>For each frame, every ATG that is not already inside an ORF in that frame starts a new ORF,
 * which ends at the first in-frame stop codon. ORFs without a stop codon are not reported.
 */
public final class OrfFinder {

    private OrfFinder() {
    }

    public static List<Orf> find(String seq, int minProteinLength) {
        String fwd = seq.toUpperCase();
        String rev = DnaTools.reverseComplement(fwd);
        int n = fwd.length();
        List<Orf> orfs = new ArrayList<>();

        for (int frame = 0; frame < 3; frame++) {
            final int f = frame;
            scan(fwd, frame, minProteinLength, (s, e, protein) ->
                    orfs.add(new Orf('+', f + 1, s + 1, e, protein)));
            scan(rev, frame, minProteinLength, (s, e, protein) ->
                    // map reverse-strand coordinates [s, e) back onto the forward strand
                    orfs.add(new Orf('-', f + 1, n - e + 1, n - s, protein)));
        }
        orfs.sort(Comparator.comparingInt(Orf::lengthNt).reversed()
                .thenComparingInt(Orf::start));
        return orfs;
    }

    private interface Sink {
        void accept(int start0, int endExclusive, String protein);
    }

    private static void scan(String s, int frame, int minLen, Sink sink) {
        int i = frame;
        while (i + 3 <= s.length()) {
            if (s.startsWith("ATG", i)) {
                StringBuilder protein = new StringBuilder();
                int j = i;
                boolean closed = false;
                while (j + 3 <= s.length()) {
                    char aa = Translator.translateCodon(s.substring(j, j + 3));
                    j += 3;
                    if (aa == '*') {
                        closed = true;
                        break;
                    }
                    protein.append(aa);
                }
                if (!closed) {
                    return; // runs off the end: no more complete ORFs in this frame
                }
                if (protein.length() >= minLen) {
                    sink.accept(i, j, protein.toString());
                }
                i = j; // continue after the stop codon
            } else {
                i += 3;
            }
        }
    }
}
