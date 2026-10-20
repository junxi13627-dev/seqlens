package io.github.seqlens;

import java.util.Map;
import java.util.TreeMap;

/**
 * Core operations on DNA sequences.
 */
public final class DnaTools {

    private DnaTools() {
    }

    /** Returns true if the sequence contains only IUPAC nucleotide codes. */
    public static boolean isValidDna(String seq) {
        return seq.toUpperCase().matches("[ACGTNRYSWKMBDHV]*");
    }

    /** Counts each base (A, C, G, T, and N/other). */
    public static Map<Character, Integer> baseCounts(String seq) {
        Map<Character, Integer> counts = new TreeMap<>();
        for (char c : "ACGTN".toCharArray()) {
            counts.put(c, 0);
        }
        for (char c : seq.toUpperCase().toCharArray()) {
            char key = "ACGT".indexOf(c) >= 0 ? c : 'N';
            counts.merge(key, 1, Integer::sum);
        }
        return counts;
    }

    /**
     * GC content as a fraction of unambiguous bases (A/C/G/T). Ambiguous bases are ignored.
     * Returns 0 for a sequence with no unambiguous bases.
     */
    public static double gcContent(String seq) {
        int gc = 0;
        int total = 0;
        for (char c : seq.toUpperCase().toCharArray()) {
            switch (c) {
                case 'G', 'C' -> { gc++; total++; }
                case 'A', 'T' -> total++;
                default -> { }
            }
        }
        return total == 0 ? 0.0 : (double) gc / total;
    }

    /** Complement of a single IUPAC nucleotide code. */
    public static char complement(char base) {
        return switch (Character.toUpperCase(base)) {
            case 'A' -> 'T';
            case 'T', 'U' -> 'A';
            case 'G' -> 'C';
            case 'C' -> 'G';
            case 'R' -> 'Y';
            case 'Y' -> 'R';
            case 'K' -> 'M';
            case 'M' -> 'K';
            case 'B' -> 'V';
            case 'V' -> 'B';
            case 'D' -> 'H';
            case 'H' -> 'D';
            case 'S' -> 'S';
            case 'W' -> 'W';
            default -> 'N';
        };
    }

    public static String reverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(complement(seq.charAt(i)));
        }
        return sb.toString();
    }

    /** DNA -> mRNA (coding strand convention: T becomes U). */
    public static String transcribe(String seq) {
        return seq.toUpperCase().replace('T', 'U');
    }

    /** Approximate melting temperature (°C). Wallace rule for &lt; 14 nt, GC-based formula otherwise. */
    public static double meltingTemperature(String seq) {
        Map<Character, Integer> c = baseCounts(seq);
        int a = c.get('A'), t = c.get('T'), g = c.get('G'), cc = c.get('C');
        int n = a + t + g + cc;
        if (n == 0) {
            return 0.0;
        }
        if (n < 14) {
            return 2.0 * (a + t) + 4.0 * (g + cc);
        }
        return 64.9 + 41.0 * (g + cc - 16.4) / n;
    }

    /** Counts all overlapping k-mers made only of A/C/G/T. */
    public static Map<String, Integer> kmerCounts(String seq, int k) {
        if (k <= 0) {
            throw new IllegalArgumentException("k must be positive");
        }
        String s = seq.toUpperCase();
        Map<String, Integer> counts = new TreeMap<>();
        for (int i = 0; i + k <= s.length(); i++) {
            String kmer = s.substring(i, i + k);
            if (kmer.matches("[ACGT]+")) {
                counts.merge(kmer, 1, Integer::sum);
            }
        }
        return counts;
    }
}
