package io.github.seqlens;

/**
 * A single entry from a FASTA file: a header line (without the leading '>') and its sequence.
 */
public record FastaRecord(String id, String description, String sequence) {

    public FastaRecord {
        if (id == null || id.isBlank()) {
            throw new IllegalArgumentException("FASTA record id must not be blank");
        }
        sequence = sequence == null ? "" : sequence.toUpperCase();
        description = description == null ? "" : description;
    }

    public int length() {
        return sequence.length();
    }
}
