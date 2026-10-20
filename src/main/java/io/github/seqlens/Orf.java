package io.github.seqlens;

/**
 * An open reading frame.
 *
 * @param strand  '+' for the forward strand, '-' for the reverse complement
 * @param frame   reading frame 1-3 on its strand
 * @param start   1-based start position on the forward strand (inclusive)
 * @param end     1-based end position on the forward strand (inclusive, includes the stop codon)
 * @param protein translated protein, without the trailing stop
 */
public record Orf(char strand, int frame, int start, int end, String protein) {

    public int lengthNt() {
        return end - start + 1;
    }
}
