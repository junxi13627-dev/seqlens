package io.github.seqlens;

import java.util.HashMap;
import java.util.Map;

/**
 * Translates DNA into protein using the standard genetic code (NCBI table 1).
 */
public final class Translator {

    private static final String BASES = "TCAG";
    // Amino acids in TCAG x TCAG x TCAG order (NCBI standard code). '*' = stop.
    private static final String AMINO =
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    private static final Map<String, Character> CODON_TABLE = new HashMap<>();

    static {
        int i = 0;
        for (char b1 : BASES.toCharArray()) {
            for (char b2 : BASES.toCharArray()) {
                for (char b3 : BASES.toCharArray()) {
                    CODON_TABLE.put("" + b1 + b2 + b3, AMINO.charAt(i++));
                }
            }
        }
    }

    private Translator() {
    }

    /** Translates a single codon; unknown or ambiguous codons give 'X'. Accepts U for T. */
    public static char translateCodon(String codon) {
        String c = codon.toUpperCase().replace('U', 'T');
        return CODON_TABLE.getOrDefault(c, 'X');
    }

    /**
     * Translates a sequence in the given reading frame (0, 1 or 2), reading through stop codons
     * (shown as '*'). Trailing bases that don't form a full codon are ignored.
     */
    public static String translate(String seq, int frame) {
        if (frame < 0 || frame > 2) {
            throw new IllegalArgumentException("frame must be 0, 1 or 2");
        }
        StringBuilder protein = new StringBuilder();
        for (int i = frame; i + 3 <= seq.length(); i += 3) {
            protein.append(translateCodon(seq.substring(i, i + 3)));
        }
        return protein.toString();
    }

    public static String translate(String seq) {
        return translate(seq, 0);
    }
}
