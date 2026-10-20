package io.github.seqlens;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Parses FASTA-formatted text into {@link FastaRecord}s.
 *
 * <p>Supports multi-line sequences, blank lines, and ';' comment lines.
 */
public final class FastaReader {

    private FastaReader() {
    }

    public static List<FastaRecord> read(Path path) throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(reader);
        }
    }

    public static List<FastaRecord> parse(String text) {
        try {
            return read(new StringReader(text));
        } catch (IOException e) {
            throw new IllegalStateException("Unexpected I/O error reading a string", e);
        }
    }

    public static List<FastaRecord> read(Reader source) throws IOException {
        BufferedReader reader = source instanceof BufferedReader br ? br : new BufferedReader(source);
        List<FastaRecord> records = new ArrayList<>();
        String header = null;
        StringBuilder seq = new StringBuilder();
        String line;
        int lineNo = 0;

        while ((line = reader.readLine()) != null) {
            lineNo++;
            line = line.strip();
            if (line.isEmpty() || line.startsWith(";")) {
                continue;
            }
            if (line.startsWith(">")) {
                if (header != null) {
                    records.add(build(header, seq));
                }
                header = line.substring(1).strip();
                seq.setLength(0);
            } else {
                if (header == null) {
                    throw new IOException("Line " + lineNo + ": sequence data found before any '>' header");
                }
                for (char c : line.toCharArray()) {
                    if (!Character.isWhitespace(c)) {
                        seq.append(c);
                    }
                }
            }
        }
        if (header != null) {
            records.add(build(header, seq));
        }
        return records;
    }

    private static FastaRecord build(String header, StringBuilder seq) {
        String[] parts = header.split("\\s+", 2);
        String id = parts[0].isEmpty() ? "unnamed" : parts[0];
        String desc = parts.length > 1 ? parts[1] : "";
        return new FastaRecord(id, desc, seq.toString());
    }
}
