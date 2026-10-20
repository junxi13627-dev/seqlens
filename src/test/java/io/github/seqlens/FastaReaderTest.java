package io.github.seqlens;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.io.StringReader;
import java.util.List;
import org.junit.jupiter.api.Test;

class FastaReaderTest {

    @Test
    void parsesMultipleMultiLineRecords() {
        String text = """
                >seq1 first test sequence
                ACGT
                acgt

                ; a comment line
                >seq2
                GGGG
                """;
        List<FastaRecord> records = FastaReader.parse(text);

        assertEquals(2, records.size());
        assertEquals("seq1", records.get(0).id());
        assertEquals("first test sequence", records.get(0).description());
        assertEquals("ACGTACGT", records.get(0).sequence());
        assertEquals("seq2", records.get(1).id());
        assertEquals("GGGG", records.get(1).sequence());
    }

    @Test
    void emptyInputGivesNoRecords() {
        assertEquals(0, FastaReader.parse("").size());
    }

    @Test
    void sequenceBeforeHeaderIsAnError() {
        assertThrows(IOException.class, () -> FastaReader.read(new StringReader("ACGT\n>x\nA")));
    }
}
